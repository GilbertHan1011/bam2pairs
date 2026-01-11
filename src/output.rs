use crate::flash::PairOutput;
use crate::segment::Statistics;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

/// Write pairs to output file
pub fn write_pairs<P: AsRef<Path>>(
    pairs: &[PairOutput],
    output_path: P,
    read_coords: bool,
    extra_tag_name: Option<&str>,
) -> Result<(), std::io::Error> {
    let file = File::create(output_path)?;
    let mut writer = BufWriter::new(file);
    
    // Write header
    let mut header = "#readID\tchr1\tpos1\tchr2\tpos2\tstrand1\tstrand2\tpair_type\tmapq1\tmapq2".to_string();
    
    if read_coords {
        header.push_str("\tRG\tstart1\tend1\tstart2\tend2");
    }
    
    if let Some(tag_name) = extra_tag_name {
        header.push_str("\t");
        header.push_str(tag_name);
    }
    
    writeln!(writer, "{}", header)?;
    
    // Write pairs
    for pair in pairs {
        writeln!(writer, "{}", pair.to_string(read_coords, extra_tag_name))?;
    }
    
    writer.flush()?;
    Ok(())
}

/// Write statistics log file
pub fn write_log<P: AsRef<Path>>(
    stats: &Statistics,
    output_path: P,
) -> Result<(), std::io::Error> {
    let file = File::create(output_path)?;
    let mut writer = BufWriter::new(file);
    
    writeln!(writer, "lowMap\t{}", stats.low_map)?;
    writeln!(writer, "manyHits\t{}", stats.many_hits)?;
    writeln!(writer, "unpaired\t{}", stats.unpaired)?;
    writeln!(writer, "selfCircle\t{}", stats.self_circle)?;
    writeln!(writer, "trans\t{}", stats.trans)?;
    writeln!(writer, "cis10K\t{}", stats.cis10k)?;
    writeln!(writer, "cis1K\t{}", stats.cis1k)?;
    writeln!(writer, "cis0\t{}", stats.cis0)?;
    writeln!(writer, "avgMapQ\t{:.2}", stats.avg_mapq())?;
    writeln!(writer, "pairTypeUU\t{}", stats.pair_type_uu)?;
    writeln!(writer, "pairTypeNU\t{}", stats.pair_type_nu)?;
    writeln!(writer, "pairTypeNN\t{}", stats.pair_type_nn)?;
    writeln!(writer, "pairTypeMM\t{}", stats.pair_type_mm)?;
    writeln!(writer, "pairTypeUR\t{}", stats.pair_type_ur)?;
    writeln!(writer, "pairTypeSC\t{}", stats.pair_type_sc)?;
    
    writer.flush()?;
    Ok(())
}

/// Write SAM records to output file
pub fn write_sam<P: AsRef<Path>>(
    records: &[Vec<u8>],
    output_path: P,
) -> Result<(), std::io::Error> {
    let file = File::create(output_path)?;
    let mut writer = BufWriter::new(file);
    
    for record in records {
        writer.write_all(record)?;
        writeln!(writer)?;
    }
    
    writer.flush()?;
    Ok(())
}
