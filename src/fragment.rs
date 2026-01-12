use bed_utils::bed::{io::Reader, BEDLike, BED};
use bed_utils::intervaltree::{Interval, Lapper};
use std::collections::HashMap;
use std::fs::File;
use std::path::Path;

/// Represents a restriction fragment
#[derive(Debug, Clone)]
pub struct Fragment {
    pub chr: String,
    pub start: u64,
    pub end: u64,
    pub name: String,
}

/// Map of restriction fragments organized by chromosome using interval trees
#[derive(Debug)]
pub struct FragmentMap {
    lapper_map: HashMap<String, Lapper<u64, BED<6>>>,
}

impl FragmentMap {
    /// Load restriction fragments from a BED file
    pub fn load<P: AsRef<Path>>(
        path: P,
        min_frag_size: Option<u64>,
        max_frag_size: Option<u64>,
    ) -> Result<Self, Box<dyn std::error::Error>> {
        let bed_file = File::open(path)?;
        let bed_reader = Reader::new(bed_file, None);
        let bed_records: Vec<BED<6>> = bed_reader
            .into_records::<BED<6>>()
            .map(|r| r.unwrap())
            .collect();

        // Filter fragments by size if specified
        let filtered_records: Vec<BED<6>> = bed_records
            .into_iter()
            .filter(|bed| {
                let frag_len = bed.end() - bed.start();
                if let Some(min_size) = min_frag_size {
                    if frag_len < min_size {
                        return false;
                    }
                }
                if let Some(max_size) = max_frag_size {
                    if frag_len > max_size {
                        return false;
                    }
                }
                true
            })
            .collect();

        let mut chrom_to_interval: HashMap<String, Vec<Interval<u64, BED<6>>>> = HashMap::new();

        for bed in &filtered_records {
            let chrom = bed.chrom().to_string();
            let interval = Interval {
                start: bed.start(),
                stop: bed.end(),
                val: bed.clone(),
            };

            chrom_to_interval
                .entry(chrom)
                .or_insert_with(Vec::new)
                .push(interval);
        }

        let lapper_map = chrom_to_interval
            .into_iter()
            .map(|(chrom, intervals)| (chrom, Lapper::new(intervals)))
            .collect();

        Ok(FragmentMap { lapper_map })
    }

    /// Get chromosome name variants to try (with and without "chr" prefix)
    fn get_chrom_variants(chrom: &str) -> Vec<String> {
        let mut variants = vec![chrom.to_string()];
        if chrom.starts_with("chr") {
            // Try without "chr" prefix
            if let Some(stripped) = chrom.strip_prefix("chr") {
                variants.push(stripped.to_string());
            }
        } else {
            // Try with "chr" prefix
            variants.push(format!("chr{}", chrom));
        }
        variants
    }

    /// Query the fragment that overlaps with a given position
    /// Uses the middle position of the read alignment span (matching HiC-Pro behavior)
    /// Tries both chromosome name formats (with and without "chr" prefix) to handle mismatches
    pub fn query_fragment(&self, chrom: &str, pos: usize) -> Option<Fragment> {
        let pos_u64 = pos as u64;
        
        // Try different chromosome name formats to handle "chr2" vs "2" mismatches
        for chrom_variant in Self::get_chrom_variants(chrom) {
            if let Some(lapper) = self.lapper_map.get(&chrom_variant) {
                let overlapping_frag: Vec<_> = lapper.find(pos_u64, pos_u64 + 1).collect();
                if overlapping_frag.len() == 1 {
                    let bed = &overlapping_frag[0].val;
                    return Some(Fragment {
                        chr: bed.chrom().to_string(),
                        start: bed.start(),
                        end: bed.end(),
                        name: bed.name().unwrap_or("").to_string(),
                    });
                } else if overlapping_frag.len() > 1 {
                    // Multiple fragments - ambiguous, return None
                    return None;
                }
            }
        }
        None
    }

    /// Query fragment using a range (start, end) instead of single position
    /// Tries both chromosome name formats (with and without "chr" prefix)
    pub fn query_fragment_range(&self, chrom: &str, start: usize, end: usize) -> Option<Fragment> {
        let start_u64 = start as u64;
        let end_u64 = end as u64;
        
        // Try different chromosome name formats
        for chrom_variant in Self::get_chrom_variants(chrom) {
            if let Some(lapper) = self.lapper_map.get(&chrom_variant) {
                let overlapping_frag: Vec<_> = lapper.find(start_u64, end_u64).collect();
                if overlapping_frag.len() == 1 {
                    let bed = &overlapping_frag[0].val;
                    return Some(Fragment {
                        chr: bed.chrom().to_string(),
                        start: bed.start(),
                        end: bed.end(),
                        name: bed.name().unwrap_or("").to_string(),
                    });
                } else if overlapping_frag.len() > 1 {
                    // Multiple fragments - ambiguous, return None
                    return None;
                }
            }
        }
        None
    }
}

/// Helper to convert BED<6> to Fragment
impl From<&BED<6>> for Fragment {
    fn from(bed: &BED<6>) -> Self {
        Fragment {
            chr: bed.chrom().to_string(),
            start: bed.start(),
            end: bed.end(),
            name: bed.name().unwrap_or("").to_string(),
        }
    }
}
