mod cigar;
mod config;
mod flash;
mod fragment;
mod grouper;
mod hic_classifier;
mod integrity;
mod output;
mod segment;
mod unc;

use clap::Parser;
use config::Config;
use flash::PairOutput;
use grouper::ReadGrouper;
use rayon::prelude::*;
use segment::Statistics;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::path::PathBuf;
use std::sync::{Arc, Mutex};
use noodles::sam::alignment::RecordBuf;

#[derive(Parser, Debug)]
#[command(name = "bam2pairs")]
#[command(about = "Extract Hi-C contact pairs from BAM files", long_about = None)]
struct Args {
    /// Input BAM file
    #[arg(short, long)]
    input: PathBuf,

    /// Processing mode: flash or unc
    #[arg(short, long, value_parser = ["flash", "unc"])]
    mode: String,

    /// Output prefix
    #[arg(short, long)]
    output: PathBuf,

    /// Number of threads
    #[arg(short, long, default_value_t = 4)]
    threads: usize,

    /// Minimum mapped ratio
    #[arg(long, default_value_t = 0.5)]
    min_mapped_ratio: f32,

    /// Minimum mapping quality
    #[arg(long, default_value_t = 10)]
    min_mapq: u8,

    /// Write SAM output
    #[arg(long)]
    write_sam: bool,

    /// Output read start and end coordinates
    #[arg(long)]
    read_coords: bool,

    /// Extract specific BAM tag (e.g., RG) to output
    #[arg(long)]
    tag: Option<String>,

    /// Restriction fragment BED file for HiC interaction classification
    #[arg(long)]
    fragment_bed: Option<PathBuf>,

    /// Minimum fragment size to consider
    #[arg(long)]
    min_frag_size: Option<u64>,

    /// Maximum fragment size to consider
    #[arg(long)]
    max_frag_size: Option<u64>,

    /// Minimum distance for cis interactions
    #[arg(long)]
    min_cis_dist: Option<u64>,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();

    // Set up Rayon thread pool
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()?;

    // Load restriction fragments if provided
    let fragment_map = if let Some(ref bed_path) = args.fragment_bed {
        eprintln!("INFO: Loading restriction fragments from: {}", bed_path.display());
        let map = fragment::FragmentMap::load(
            bed_path,
            args.min_frag_size,
            args.max_frag_size,
        )?;
        eprintln!("INFO: Loaded restriction fragments");
        Some(Arc::new(map))
    } else {
        None
    };

    // Open BAM file
    let file = File::open(&args.input)?;
    let reader = BufReader::new(file);
    let mut bam_reader = noodles::bam::io::Reader::new(reader);

    // Read header
    let sam_header = bam_reader.read_header()?;

    // Check BAM sort order: bam2pairs requires queryname-sorted input.
    //
    // Note: Many aligners output coordinate-sorted BAM by default; grouping by read name
    // (and thus correct pairing / multi-hit handling) requires name-sorted BAM.
    //
    // The noodles SAM header API in this version does not expose a direct `sort_order`
    // accessor on the top-level header map, so we conservatively inspect the debug
    // representation of the header for an SO:coordinate flag.
    let header_text = format!("{sam_header:?}");
    if header_text.contains("SO:Coordinate") {
        eprintln!("\x1b[31mERROR: Input BAM file is coordinate-sorted!\x1b[0m");
        eprintln!("       bam2pairs requires a name-sorted BAM (QueryName) to correctly group paired-end reads.");
        eprintln!("       Please run: \x1b[1msamtools sort -n -o <out_name_sorted.bam> <in.bam>\x1b[0m");
        return Err("Input BAM is coordinate-sorted. Aborting.".into());
    }

    // Build reference sequence ID to chromosome name mapping
    let mut ref_seq_map = HashMap::new();
    for (idx, (name, _)) in sam_header.reference_sequences().iter().enumerate() {
        if let Ok(name_str) = std::str::from_utf8(name.as_ref()) {
            ref_seq_map.insert(idx, name_str.to_string());
        }
    }
    let ref_seq_map = Arc::new(ref_seq_map);

    // Create configuration
    let config = Config::new(
        args.min_mapq,
        args.min_mapped_ratio,
        args.tag.clone(),
        fragment_map,
        args.min_cis_dist,
        ref_seq_map,
    );

    eprintln!("INFO: Reading BAM file: {}", args.input.display());
    eprintln!("INFO: Mode: {}", args.mode);
    eprintln!("INFO: Threads: {}", args.threads);
    eprintln!("INFO: Min MapQ: {}", config.min_mapq);
    eprintln!("INFO: Min mapped ratio: {}", config.min_mapped_ratio);
    if args.read_coords {
        eprintln!("INFO: Output read coordinates: enabled");
    }
    if let Some(tag) = &config.extract_tag {
        eprintln!("INFO: Extracting Tag: {}", tag);
    }
    if config.fragment_map.is_some() {
        eprintln!("INFO: HiC interaction classification: enabled");
        if let Some(min_dist) = config.min_cis_dist {
            eprintln!("INFO: Minimum cis distance: {}", min_dist);
        }
    }

    eprintln!("INFO: Grouping reads by name (streaming)...");
    let mut grouper = ReadGrouper::new(bam_reader, sam_header);

    // Prepare output writer with header
    let pairs_file = format!("{}.{}.pairs", args.output.display(), args.mode);
    eprintln!("INFO: Writing pairs to {}", pairs_file);
    let pairs_fh = File::create(&pairs_file)?;
    let mut pairs_writer = BufWriter::new(pairs_fh);

    let mut header_line = "#readID\tchr1\tpos1\tchr2\tpos2\tstrand1\tstrand2\tpair_type\tmapq1\tmapq2".to_string();
    if args.read_coords {
        header_line.push_str("\tstart1\tend1\tstart2\tend2");
    }
    if let Some(tag_name) = config.extract_tag.as_deref() {
        header_line.push_str("\t");
        header_line.push_str(tag_name);
    }
    if config.fragment_map.is_some() {
        header_line.push_str("\thic_type");
    }
    writeln!(pairs_writer, "{}", header_line)?;

    // Streaming processing with chunked parallelism
    eprintln!("INFO: Processing pairs in streaming mode...");
    let mut stats = Statistics::new();
    let chunk_size: usize = 10_000;
    let mut chunk: Vec<Vec<RecordBuf>> = Vec::with_capacity(chunk_size);
    let include_hic_type = config.fragment_map.is_some();

    while let Some(group) = grouper.next_group() {
        chunk.push(group);

        if chunk.len() >= chunk_size {
            let (pairs, chunk_stats) =
                process_chunk_parallel(&chunk, &config, args.mode.as_str(), include_hic_type);
            stats.merge(&chunk_stats);

            for pair in pairs {
                writeln!(
                    pairs_writer,
                    "{}",
                    pair.to_string(args.read_coords, config.extract_tag.as_deref(), include_hic_type)
                )?;
            }

            chunk.clear();
        }
    }

    // Process any remaining groups
    if !chunk.is_empty() {
        let (pairs, chunk_stats) =
            process_chunk_parallel(&chunk, &config, args.mode.as_str(), include_hic_type);
        stats.merge(&chunk_stats);

        for pair in pairs {
            writeln!(
                pairs_writer,
                "{}",
                pair.to_string(args.read_coords, config.extract_tag.as_deref(), include_hic_type)
            )?;
        }
    }

    pairs_writer.flush()?;

    eprintln!("INFO: Processed {} valid pairs", stats.valid_pairs);

    let log_file = format!("{}.{}2pairs.log", args.output.display(), args.mode);
    eprintln!("INFO: Writing statistics to {}", log_file);
    output::write_log(&stats, &log_file)?;

    // Write SAM if requested
    if args.write_sam {
        let sam_file = format!("{}.{}.sam", args.output.display(), args.mode);
        eprintln!("INFO: Writing SAM to {}", sam_file);
        // SAM writing would require storing original records
        eprintln!("WARN: SAM output not yet implemented");
    }

    eprintln!("INFO: Done!");
    Ok(())
}

fn process_chunk_parallel(
    chunk: &[Vec<RecordBuf>],
    config: &Config,
    mode: &str,
    include_hic_type: bool,
) -> (Vec<PairOutput>, Statistics) {
    let pairs_mutex: Mutex<Vec<PairOutput>> = Mutex::new(Vec::new());
    let stats_mutex: Mutex<Statistics> = Mutex::new(Statistics::new());

    chunk.par_iter().for_each(|group| {
        let mut local_stats = Statistics::new();

        let pair_result = match mode {
            "flash" => crate::flash::process_flash(group, config, &mut local_stats),
            "unc" => crate::unc::process_unc(group, config, &mut local_stats),
            _ => None,
        };

        if let Some(pair) = pair_result {
            // Only keep pairs; HiC type string is already embedded if enabled
            pairs_mutex.lock().unwrap().push(pair);
        }

        stats_mutex.lock().unwrap().merge(&local_stats);
    });

    let pairs = pairs_mutex.into_inner().unwrap();
    let stats = stats_mutex.into_inner().unwrap();

    (pairs, stats)
}
