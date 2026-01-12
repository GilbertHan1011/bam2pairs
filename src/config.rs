use crate::fragment::FragmentMap;
use std::collections::HashMap;
use std::sync::Arc;

/// Configuration constants for pair processing
#[derive(Debug, Clone)]
pub struct Config {
    pub min_mapq: u8,
    pub min_mapped_ratio: f32,
    pub min_clip_size: usize,
    pub max_self_circle_dist: usize,
    pub max_pair_dist: usize,
    pub extract_tag: Option<String>,
    pub fragment_map: Option<Arc<FragmentMap>>,
    pub min_cis_dist: Option<u64>,
    pub ref_seq_map: Arc<HashMap<usize, String>>, // Map ref seq ID to chromosome name
}

impl Config {
    pub fn new(
        min_mapq: u8,
        min_mapped_ratio: f32,
        extract_tag: Option<String>,
        fragment_map: Option<Arc<FragmentMap>>,
        min_cis_dist: Option<u64>,
        ref_seq_map: Arc<HashMap<usize, String>>,
    ) -> Self {
        Self {
            min_mapq,
            min_mapped_ratio,
            min_clip_size: 20,
            max_self_circle_dist: 10,
            max_pair_dist: 1000,
            extract_tag,
            fragment_map,
            min_cis_dist,
            ref_seq_map,
        }
    }
    
    /// Get chromosome name from reference sequence ID
    pub fn get_chrom_name(&self, ref_seq_id: Option<usize>) -> String {
        ref_seq_id
            .and_then(|id| self.ref_seq_map.get(&id).cloned())
            .unwrap_or_else(|| "!".to_string())
    }
}

impl Default for Config {
    fn default() -> Self {
        Self::new(10, 0.5, None, None, None, Arc::new(HashMap::new()))
    }
}
