use crate::cigar::cigar_to_segment;
use crate::config::Config;
use crate::hic_classifier::{classify_interaction, ReadInfo as HicReadInfo};
use crate::integrity::{check_integrity_1_seg, check_integrity_2_seg};
use crate::segment::{Segment, Statistics};
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record_buf::data::field::Value;
use noodles::sam::alignment::RecordBuf;

#[derive(Debug, Clone)]
pub struct PairOutput {
    pub read_id: String,
    pub chr1: String,
    pub pos1: usize,
    pub chr2: String,
    pub pos2: usize,
    pub strand1: char,
    pub strand2: char,
    pub pair_type: String,
    pub mapq1: u8,
    pub mapq2: u8,
    // Optional extended fields
    pub rg: Option<String>,
    pub start1: usize,
    pub end1: usize,
    pub start2: usize,
    pub end2: usize,
    pub extra_tag: Option<String>,
    pub hic_type: Option<String>, // HiC interaction type (VI, DE, RE, SC, etc.)
}

impl PairOutput {
    pub fn to_string(&self, show_read_coords: bool, _extra_tag_name: Option<&str>, include_hic_type: bool) -> String {
        let basic = format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.read_id,
            self.chr1,
            self.pos1,
            self.chr2,
            self.pos2,
            self.strand1,
            self.strand2,
            self.pair_type,
            self.mapq1,
            self.mapq2
        );

        let mut result = basic;
        
        if show_read_coords {
            result = format!(
                "{}\t{}\t{}\t{}\t{}",
                result,
                self.start1,
                self.end1,
                self.start2,
                self.end2
            );
        }
        
        if let Some(tag_val) = &self.extra_tag {
            result = format!("{}\t{}", result, tag_val);
        }
        
        if include_hic_type {
            let hic_type_str = self.hic_type.as_deref().unwrap_or("*");
            result = format!("{}\t{}", result, hic_type_str);
        }
        
        result
    }
}

/// Helper to extract Read Group (RG) tag
fn get_rg(record: &RecordBuf) -> Option<String> {
    record
        .data()
        .get(&Tag::READ_GROUP)
        .and_then(|val| {
            match val {
                Value::String(s) => Some(s.to_string()),
                _ => None,
            }
        })
}

/// Helper to extract a specific tag from a record
fn get_tag_value(record: &RecordBuf, tag_str: &Option<String>) -> Option<String> {
    let tag_name = tag_str.as_ref()?;
    let bytes = tag_name.as_bytes();
    
    // Safety check for tag length
    if bytes.len() != 2 {
        return None;
    }
    
    let tag = Tag::new(bytes[0], bytes[1]);

    // .get() returns Option<&Value>
    // Match on Value enum variants to extract the actual value
    record.data().get(&tag).and_then(|val| {
        match val {
            Value::String(s) => Some(s.to_string()),
            Value::Int8(i) => Some(i.to_string()),
            Value::UInt8(u) => Some(u.to_string()),
            Value::Int16(i) => Some(i.to_string()),
            Value::UInt16(u) => Some(u.to_string()),
            Value::Int32(i) => Some(i.to_string()),
            Value::UInt32(u) => Some(u.to_string()),
            Value::Float(f) => Some(f.to_string()),
            _ => Some(".".to_string()),
        }
    })
}

/// Process flash mode (merged reads)
/// Port of flash2pairs from flash2pairs.h
pub fn process_flash(
    records: &[RecordBuf],
    config: &Config,
    stats: &mut Statistics,
) -> Option<PairOutput> {
    if records.is_empty() {
        return None;
    }
    
    if records.len() == 1 {
        // Single record: intra-read contact
        process_single_record(&records[0], config, stats)
    } else if records.len() == 2 {
        // Two records: chimeric alignment
        process_two_records(&records[0], &records[1], config, stats)
    } else {
        // Too many segments
        stats.many_hits += 1;
        None
    }
}

fn process_single_record(
    record: &RecordBuf,
    config: &Config,
    stats: &mut Statistics,
) -> Option<PairOutput> {
    let read_id = String::from_utf8_lossy(record.name().unwrap()).to_string();
    let flags = record.flags();
    let mapq = record.mapping_quality().map(|q| q.get()).unwrap_or(0);
    
    // Check if mapped
    let is_mapped = !flags.is_unmapped()
        && record.reference_sequence_id().is_some()
        && !record.cigar().as_ref().is_empty();
    
    let chr1 = if is_mapped {
        config.get_chrom_name(record.reference_sequence_id())
    } else {
        "!".to_string()
    };
    
    let pos1 = record.alignment_start().map(|p| p.get()).unwrap_or(0);
    let mut pair_type = if is_mapped { "UU" } else { "NN" };
    
    let segment = if is_mapped {
        match cigar_to_segment(record.cigar(), pos1) {
            Ok(seg) => seg,
            Err(_) => {
                stats.many_hits += 1;
                return None;
            }
        }
    } else {
        Segment::new()
    };
    
    if is_mapped && segment.seg_cnt > 2 {
        stats.many_hits += 1;
        stats.pair_type_mm += 1;
        pair_type = "MM";
    }
    
    if is_mapped && !check_integrity_1_seg(&segment, config) {
        stats.low_map += 1;
        pair_type = "NN";
    }
    
    let pos2 = if is_mapped && segment.seg_cnt > 0 {
        segment.right[segment.seg_cnt - 1]
    } else {
        pos1
    };
    
    let dist = if pos2 > pos1 { pos2 - pos1 } else { 0 };
    
    if is_mapped {
        if dist >= 10000 {
            stats.cis10k += 1;
        } else if dist >= 1000 {
            stats.cis1k += 1;
        } else {
            stats.cis0 += 1;
        }
    }
    
    stats.sum_mapq += mapq as u64;
    stats.valid_pairs += 1;
    
    match pair_type {
        "UU" => stats.pair_type_uu += 1,
        "NN" => stats.pair_type_nn += 1,
        "MM" => stats.pair_type_mm += 1,
        _ => {}
    }
    
    let rg = get_rg(record);
    let (start1, end1) = segment.span();
    // Since it's a single read, both "sides" have the same read span
    let (start2, end2) = (start1, end1);
    let extra_tag = get_tag_value(record, &config.extract_tag);
    
    // Classify HiC interaction if fragment map is available
    let hic_type = if let Some(ref frag_map) = config.fragment_map {
        // Calculate middle position for fragment query (matching HiC-Pro behavior)
        let middle_pos1 = if is_mapped && segment.seg_cnt > 0 {
            (start1 + end1) / 2
        } else {
            pos1
        };
        let middle_pos2 = if is_mapped && segment.seg_cnt > 0 {
            (start2 + end2) / 2
        } else {
            pos2
        };
        
        let frag1 = if is_mapped && chr1 != "!" {
            frag_map.query_fragment(&chr1, middle_pos1)
        } else {
            None
        };
        let frag2 = if is_mapped && chr1 != "!" {
            frag_map.query_fragment(&chr1, middle_pos2)
        } else {
            None
        };
        
        let r1_info = HicReadInfo {
            chr: chr1.clone(),
            pos: pos1,
            strand: '+',
            is_mapped,
        };
        let r2_info = HicReadInfo {
            chr: chr1.clone(),
            pos: pos2,
            strand: '-',
            is_mapped,
        };
        
        Some(classify_interaction(
            &r1_info,
            &r2_info,
            frag1.as_ref(),
            frag2.as_ref(),
            config.min_cis_dist,
        ).as_str().to_string())
    } else {
        None
    };
    
    Some(PairOutput {
        read_id,
        chr1: chr1.clone(),
        pos1,
        chr2: chr1,
        pos2,
        strand1: '+',
        strand2: '-',
        pair_type: pair_type.to_string(),
        mapq1: mapq,
        mapq2: mapq,
        rg,
        start1,
        end1,
        start2,
        end2,
        extra_tag,
        hic_type,
    })
}

fn process_two_records(
    record1: &RecordBuf,
    record2: &RecordBuf,
    config: &Config,
    stats: &mut Statistics,
) -> Option<PairOutput> {
    let read_id = String::from_utf8_lossy(record1.name().unwrap()).to_string();
    
    // Process first record
    let flags1 = record1.flags();
    let mapq1 = record1.mapping_quality().map(|q| q.get()).unwrap_or(0);
    let is_mapped1 = !flags1.is_unmapped()
        && record1.reference_sequence_id().is_some()
        && !record1.cigar().as_ref().is_empty();
    
    let chr1 = if is_mapped1 {
        config.get_chrom_name(record1.reference_sequence_id())
    } else {
        "!".to_string()
    };
    
    let mut pos1 = record1.alignment_start().map(|p| p.get()).unwrap_or(0);
    let strand1 = if flags1.is_reverse_complemented() { '-' } else { '+' };
    
    let seg1 = if is_mapped1 {
        match cigar_to_segment(record1.cigar(), pos1) {
            Ok(seg) => seg,
            Err(_) => Segment::new(),
        }
    } else {
        Segment::new()
    };
    
    // Process second record
    let flags2 = record2.flags();
    let mapq2 = record2.mapping_quality().map(|q| q.get()).unwrap_or(0);
    let is_mapped2 = !flags2.is_unmapped()
        && record2.reference_sequence_id().is_some()
        && !record2.cigar().as_ref().is_empty();
    
    let chr2 = if is_mapped2 {
        config.get_chrom_name(record2.reference_sequence_id())
    } else {
        "!".to_string()
    };
    
    let mut pos2 = record2.alignment_start().map(|p| p.get()).unwrap_or(0);
    let strand2 = if flags2.is_reverse_complemented() { '-' } else { '+' };
    
    let seg2 = if is_mapped2 {
        match cigar_to_segment(record2.cigar(), pos2) {
            Ok(seg) => seg,
            Err(_) => Segment::new(),
        }
    } else {
        Segment::new()
    };
    
    // Determine pair type
    let q1_good = is_mapped1 && mapq1 >= config.min_mapq;
    let q2_good = is_mapped2 && mapq2 >= config.min_mapq;
    
    let mut pair_type = if !q1_good && !q2_good {
        "NN"
    } else if !q1_good || !q2_good {
        "NU"
    } else {
        "UU"
    };
    
    if is_mapped1 && is_mapped2 && (seg1.seg_cnt != 1 || seg2.seg_cnt != 1) {
        stats.many_hits += 1;
        pair_type = "MM";
        stats.pair_type_mm += 1;
    }
    
    if is_mapped1 && is_mapped2 && !check_integrity_2_seg(&seg1, &seg2, config) {
        stats.low_map += 1;
        pair_type = "NN";
    }
    
    let pair_mapq = mapq1.min(mapq2);
    stats.sum_mapq += pair_mapq as u64;
    stats.valid_pairs += 1;
    
    // Get outmost positions
    if seg1.left_clip > seg1.right_clip && seg1.right.len() > 0 {
        pos1 = seg1.right[0];
    }
    if seg2.left_clip > seg2.right_clip && seg2.right.len() > 0 {
        pos2 = seg2.right[0];
    }
    
    let rg = get_rg(record1);
    let (start1, end1) = seg1.span();
    let (start2, end2) = seg2.span();
    let extra_tag = get_tag_value(record1, &config.extract_tag);
    
    // Determine order and calculate distance
    let (final_chr1, final_pos1, final_strand1, final_start1, final_end1, final_chr2, final_pos2, final_strand2, final_start2, final_end2, final_mapq1, final_mapq2) =
        if chr1 < chr2 || (chr1 == chr2 && pos1 < pos2) {
            if chr1 == chr2 {
                let dist = if pos2 > pos1 { pos2 - pos1 } else { 0 };
                if dist <= config.max_self_circle_dist {
                    stats.self_circle += 1;
                    pair_type = "SC";
                    stats.pair_type_sc += 1;
                } else {
                    if dist >= 10000 {
                        stats.cis10k += 1;
                    } else if dist >= 1000 {
                        stats.cis1k += 1;
                    } else {
                        stats.cis0 += 1;
                    }
                }
            } else {
                stats.trans += 1;
            }
            (chr1, pos1, strand1, start1, end1, chr2, pos2, strand2, start2, end2, mapq1, mapq2)
        } else {
            if chr1 == chr2 {
                let dist = if pos1 > pos2 { pos1 - pos2 } else { 0 };
                if dist <= config.max_self_circle_dist {
                    stats.self_circle += 1;
                    pair_type = "SC";
                    stats.pair_type_sc += 1;
                } else {
                    if dist >= 10000 {
                        stats.cis10k += 1;
                    } else if dist >= 1000 {
                        stats.cis1k += 1;
                    } else {
                        stats.cis0 += 1;
                    }
                }
            } else {
                stats.trans += 1;
            }
            (chr2, pos2, strand2, start2, end2, chr1, pos1, strand1, start1, end1, mapq2, mapq1)
        };
    
    // Update pair type counters
    match pair_type {
        "UU" => stats.pair_type_uu += 1,
        "NU" => stats.pair_type_nu += 1,
        "NN" => stats.pair_type_nn += 1,
        "MM" => {}, // Already counted
        "UR" => stats.pair_type_ur += 1,
        "SC" => {}, // Already counted
        _ => {}
    }
    
    // Skip self-circles
    if pair_type == "SC" {
        return None;
    }
    
    // Classify HiC interaction if fragment map is available
    let hic_type = if let Some(ref frag_map) = config.fragment_map {
        // Calculate middle position for fragment query (matching HiC-Pro behavior)
        let middle_pos1 = (final_start1 + final_end1) / 2;
        let middle_pos2 = (final_start2 + final_end2) / 2;
        
        let frag1 = if final_chr1 != "!" {
            frag_map.query_fragment(&final_chr1, middle_pos1)
        } else {
            None
        };
        let frag2 = if final_chr2 != "!" {
            frag_map.query_fragment(&final_chr2, middle_pos2)
        } else {
            None
        };
        
        let r1_info = HicReadInfo {
            chr: final_chr1.clone(),
            pos: final_pos1,
            strand: final_strand1,
            is_mapped: is_mapped1,
        };
        let r2_info = HicReadInfo {
            chr: final_chr2.clone(),
            pos: final_pos2,
            strand: final_strand2,
            is_mapped: is_mapped2,
        };
        
        Some(classify_interaction(
            &r1_info,
            &r2_info,
            frag1.as_ref(),
            frag2.as_ref(),
            config.min_cis_dist,
        ).as_str().to_string())
    } else {
        None
    };
    
    Some(PairOutput {
        read_id,
        chr1: final_chr1,
        pos1: final_pos1,
        chr2: final_chr2,
        pos2: final_pos2,
        strand1: final_strand1,
        strand2: final_strand2,
        pair_type: pair_type.to_string(),
        mapq1: final_mapq1,
        mapq2: final_mapq2,
        rg,
        start1: final_start1,
        end1: final_end1,
        start2: final_start2,
        end2: final_end2,
        extra_tag,
        hic_type,
    })
}
