/// Configuration constants for pair processing
#[derive(Debug, Clone)]
pub struct Config {
    pub min_mapq: u8,
    pub min_mapped_ratio: f32,
    pub min_clip_size: usize,
    pub max_self_circle_dist: usize,
    pub max_pair_dist: usize,
}

impl Config {
    pub fn new(min_mapq: u8, min_mapped_ratio: f32) -> Self {
        Self {
            min_mapq,
            min_mapped_ratio,
            min_clip_size: 20,
            max_self_circle_dist: 10,
            max_pair_dist: 1000,
        }
    }
}

impl Default for Config {
    fn default() -> Self {
        Self::new(10, 0.5)
    }
}
