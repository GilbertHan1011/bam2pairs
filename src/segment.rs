/// Represents a parsed segment from CIGAR string
#[derive(Debug, Clone)]
pub struct Segment {
    pub seg_cnt: usize,
    pub left_clip: usize,
    pub right_clip: usize,
    pub left: Vec<usize>,
    pub right: Vec<usize>,
    pub mappable: usize,
}

impl Segment {
    pub fn new() -> Self {
        Self {
            seg_cnt: 0,
            left_clip: 0,
            right_clip: 0,
            left: Vec::new(),
            right: Vec::new(),
            mappable: 0,
        }
    }

    /// Get the genomic span (start, end) of the segment
    pub fn span(&self) -> (usize, usize) {
        if self.seg_cnt == 0 || self.left.is_empty() {
            (0, 0)
        } else {
            (self.left[0], self.right[self.seg_cnt - 1])
        }
    }
}

impl Default for Segment {
    fn default() -> Self {
        Self::new()
    }
}

/// Represents a segment key with alignment information
#[derive(Debug, Clone)]
pub struct SegKey {
    pub chr: String,
    pub pos: usize,
    pub strand: char,
    pub cigar: String,
    pub mapq: u8,
    pub is_mapped: bool,
}

/// Statistics tracking for pair processing
#[derive(Debug, Clone, Default)]
pub struct Statistics {
    pub low_map: u64,
    pub many_hits: u64,
    pub unpaired: u64,
    pub trans: u64,
    pub self_circle: u64,
    pub cis0: u64,
    pub cis1k: u64,
    pub cis10k: u64,
    pub sum_mapq: u64,
    pub valid_pairs: u64,
    pub pair_type_uu: u64,
    pub pair_type_nu: u64,
    pub pair_type_nn: u64,
    pub pair_type_mm: u64,
    pub pair_type_ur: u64,
    pub pair_type_sc: u64,
    // HiC interaction type counters
    pub hic_valid: u64,
    pub hic_dangling_end: u64,
    pub hic_religation: u64,
    pub hic_self_circle_frag: u64,
    pub hic_single: u64,
    pub hic_filtered: u64,
    pub hic_dumped: u64,
}

impl Statistics {
    pub fn new() -> Self {
        Self::default()
    }

    /// Merge another Statistics instance into this one
    pub fn merge(&mut self, other: &Statistics) {
        self.low_map += other.low_map;
        self.many_hits += other.many_hits;
        self.unpaired += other.unpaired;
        self.trans += other.trans;
        self.self_circle += other.self_circle;
        self.cis0 += other.cis0;
        self.cis1k += other.cis1k;
        self.cis10k += other.cis10k;
        self.sum_mapq += other.sum_mapq;
        self.valid_pairs += other.valid_pairs;
        self.pair_type_uu += other.pair_type_uu;
        self.pair_type_nu += other.pair_type_nu;
        self.pair_type_nn += other.pair_type_nn;
        self.pair_type_mm += other.pair_type_mm;
        self.pair_type_ur += other.pair_type_ur;
        self.pair_type_sc += other.pair_type_sc;
        // HiC counters
        self.hic_valid += other.hic_valid;
        self.hic_dangling_end += other.hic_dangling_end;
        self.hic_religation += other.hic_religation;
        self.hic_self_circle_frag += other.hic_self_circle_frag;
        self.hic_single += other.hic_single;
        self.hic_filtered += other.hic_filtered;
        self.hic_dumped += other.hic_dumped;
    }

    /// Calculate average MapQ
    pub fn avg_mapq(&self) -> f64 {
        if self.valid_pairs > 0 {
            self.sum_mapq as f64 / self.valid_pairs as f64
        } else {
            0.0
        }
    }
}
