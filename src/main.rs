mod cigar;
mod config;
mod flash;
mod grouper;
mod integrity;
mod output;
mod segment;
mod unc;

use clap::Parser;
use config::Config;
use grouper::collect_groups;
use rayon::prelude::*;
use segment::Statistics;
use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;
use std::sync::Mutex;

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
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();

    // Set up Rayon thread pool
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()?;

    // Create configuration
    let config = Config::new(args.min_mapq, args.min_mapped_ratio);

    eprintln!("INFO: Reading BAM file: {}", args.input.display());
    eprintln!("INFO: Mode: {}", args.mode);
    eprintln!("INFO: Threads: {}", args.threads);
    eprintln!("INFO: Min MapQ: {}", config.min_mapq);
    eprintln!("INFO: Min mapped ratio: {}", config.min_mapped_ratio);

    // Open BAM file
    let file = File::open(&args.input)?;
    let reader = BufReader::new(file);
    let mut bam_reader = noodles::bam::io::Reader::new(reader);

    // Read header
    let header = bam_reader.read_header()?;

    eprintln!("INFO: Grouping reads by name...");
    let groups = collect_groups(bam_reader, header);
    eprintln!("INFO: Found {} read groups", groups.len());

    // Process groups in parallel
    eprintln!("INFO: Processing pairs...");
    let pairs_mutex = Mutex::new(Vec::new());
    let stats_mutex = Mutex::new(Statistics::new());

    groups.par_iter().for_each(|group| {
        let mut local_stats = Statistics::new();
        
        let pair_result = match args.mode.as_str() {
            "flash" => flash::process_flash(group, &config, &mut local_stats),
            "unc" => unc::process_unc(group, &config, &mut local_stats),
            _ => None,
        };

        // Collect pair if generated
        if let Some(pair) = pair_result {
            pairs_mutex.lock().unwrap().push(pair);
        }

        // Merge statistics
        stats_mutex.lock().unwrap().merge(&local_stats);
    });

    // Extract results
    let pairs = pairs_mutex.into_inner().unwrap();
    let stats = stats_mutex.into_inner().unwrap();

    eprintln!("INFO: Processed {} valid pairs", stats.valid_pairs);

    // Write output files
    let pairs_file = format!("{}.{}.pairs", args.output.display(), args.mode);
    eprintln!("INFO: Writing pairs to {}", pairs_file);
    output::write_pairs(&pairs, &pairs_file)?;

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
