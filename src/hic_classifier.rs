use crate::fragment::Fragment;

/// Hi-C interaction types
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HicType {
    VI,   // Valid Interaction
    DE,   // Dangling End
    RE,   // Religation
    SC,   // Self Circle
    SI,   // Single (one read unmapped)
    FILT, // Filtered (by distance/size criteria)
    DUMP, // Other/Unknown
}

impl HicType {
    pub fn as_str(&self) -> &'static str {
        match self {
            HicType::VI => "VI",
            HicType::DE => "DE",
            HicType::RE => "RE",
            HicType::SC => "SC",
            HicType::SI => "SI",
            HicType::FILT => "FILT",
            HicType::DUMP => "DUMP",
        }
    }
}

/// Information about a read for classification
#[derive(Debug, Clone)]
pub struct ReadInfo {
    pub chr: String,
    pub pos: usize,
    pub strand: char,
    pub is_mapped: bool,
}

/// Classify Hi-C interaction type based on reads and fragments
pub fn classify_interaction(
    r1_info: &ReadInfo,
    r2_info: &ReadInfo,
    frag1: Option<&Fragment>,
    frag2: Option<&Fragment>,
    min_cis_dist: Option<u64>,
) -> HicType {
    // Check if either read is unmapped
    if !r1_info.is_mapped || !r2_info.is_mapped {
        return HicType::SI;
    }

    // Both reads must be mapped and have fragments
    match (frag1, frag2) {
        (Some(f1), Some(f2)) => {
            if are_same_fragment(f1, f2) {
                // Same fragment - check for dangling end or self circle
                if is_self_circle(r1_info, r2_info) {
                    HicType::SC
                } else if is_dangling_end(r1_info, r2_info) {
                    HicType::DE
                } else {
                    HicType::DUMP
                }
            } else {
                // Different fragments - check for religation
                if is_religation(f1, f2) {
                    HicType::RE
                } else {
                    // Valid interaction - check distance filter
                    let cis_dist = calculate_cis_distance(r1_info, r2_info);
                    if let Some(min_dist) = min_cis_dist {
                        if let Some(dist) = cis_dist {
                            if dist < min_dist {
                                return HicType::FILT;
                            }
                        }
                    }
                    HicType::VI
                }
            }
        }
        _ => HicType::DUMP,
    }
}

/// Check if two fragments are the same
fn are_same_fragment(frag1: &Fragment, frag2: &Fragment) -> bool {
    frag1.chr == frag2.chr
        && frag1.start == frag2.start
        && frag1.end == frag2.end
}

/// Check if two fragments are contiguous (one's end == other's start)
fn are_contiguous_fragments(frag1: &Fragment, frag2: &Fragment) -> bool {
    if frag1.chr != frag2.chr {
        return false;
    }
    frag1.end == frag2.start || frag1.start == frag2.end
}

/// Check if this is a self circle interaction
/// Self circle: same fragment, ordered reads have -/+ strands
fn is_self_circle(r1_info: &ReadInfo, r2_info: &ReadInfo) -> bool {
    if r1_info.chr != r2_info.chr {
        return false;
    }
    
    // Order reads by position
    let (r1, r2) = if r1_info.pos <= r2_info.pos {
        (r1_info, r2_info)
    } else {
        (r2_info, r1_info)
    };
    
    r1.strand == '-' && r2.strand == '+'
}

/// Check if this is a dangling end interaction
/// Dangling end: same fragment, ordered reads have +/- strands
fn is_dangling_end(r1_info: &ReadInfo, r2_info: &ReadInfo) -> bool {
    if r1_info.chr != r2_info.chr {
        return false;
    }
    
    // Order reads by position
    let (r1, r2) = if r1_info.pos <= r2_info.pos {
        (r1_info, r2_info)
    } else {
        (r2_info, r1_info)
    };
    
    r1.strand == '+' && r2.strand == '-'
}

/// Check if this is a religation (contiguous fragments)
fn is_religation(frag1: &Fragment, frag2: &Fragment) -> bool {
    are_contiguous_fragments(frag1, frag2)
}

/// Calculate cis distance between two reads
fn calculate_cis_distance(r1_info: &ReadInfo, r2_info: &ReadInfo) -> Option<u64> {
    if r1_info.chr != r2_info.chr {
        return None;
    }
    
    if r1_info.pos > r2_info.pos {
        Some((r1_info.pos - r2_info.pos) as u64)
    } else {
        Some((r2_info.pos - r1_info.pos) as u64)
    }
}

/// Calculate fragment distance for different interaction types
pub fn calculate_fragment_distance(
    r1_info: &ReadInfo,
    r2_info: &ReadInfo,
    frag1: &Fragment,
    frag2: &Fragment,
    interaction_type: HicType,
) -> Option<u64> {
    // Order reads by position
    let (r1, r2, rf1, rf2) = if r1_info.pos <= r2_info.pos {
        (r1_info, r2_info, frag1, frag2)
    } else {
        (r2_info, r1_info, frag2, frag1)
    };

    match interaction_type {
        HicType::DE | HicType::RE => {
            // Distance is between read positions
            Some((r2.pos.saturating_sub(r1.pos)) as u64)
        }
        HicType::SC => {
            // Self circle: distance from read1 to fragment start + fragment end to read2
            let d1 = r1.pos.saturating_sub(rf1.start as usize);
            let d2 = (rf2.end as usize).saturating_sub(r2.pos);
            Some((d1 + d2) as u64)
        }
        HicType::VI => {
            // Valid interaction: distance from read to fragment end/start
            let dr1 = if r1.strand == '+' {
                (rf1.end as usize).saturating_sub(r1.pos)
            } else {
                r1.pos.saturating_sub(rf1.start as usize)
            };
            
            let dr2 = if r2.strand == '+' {
                (rf2.end as usize).saturating_sub(r2.pos)
            } else {
                r2.pos.saturating_sub(rf2.start as usize)
            };
            
            Some((dr1 + dr2) as u64)
        }
        _ => None,
    }
}
