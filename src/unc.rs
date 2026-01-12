use crate::cigar::cigar_to_segment;
use crate::config::Config;
use crate::flash::PairOutput;
use crate::hic_classifier::{classify_interaction, ReadInfo as HicReadInfo};
use crate::integrity::{check_integrity_1_seg, check_integrity_2_seg};
use crate::segment::{Segment, SegKey, Statistics};
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record_buf::data::field::Value;
use noodles::sam::alignment::RecordBuf;

/// Process unc mode (unmerged paired-end reads)
/// Port of unc2pairs from unc2pairs.h
pub fn process_unc(
    records: &[RecordBuf],
    config: &Config,
    stats: &mut Statistics,
) -> Option<PairOutput> {
    if records.is_empty() {
        return None;
    }
    
    // Separate R1 and R2 records
    let mut r1_records: Vec<&RecordBuf> = Vec::new();
    let mut r2_records: Vec<&RecordBuf> = Vec::new();
    
    let read_id = String::from_utf8_lossy(records[0].name().unwrap()).to_string();
    
    for record in records {
        let flags = record.flags();
        if flags.is_first_segment() {
            r1_records.push(record);
        } else if flags.is_last_segment() {
            r2_records.push(record);
        }
    }
    
    // Extract RG tag from first R1 record, or first R2 if R1 is empty
    let rg = if !r1_records.is_empty() {
        r1_records[0].data().get(&Tag::READ_GROUP).and_then(|v| {
            match v {
                Value::String(s) => Some(s.to_string()),
                _ => None,
            }
        })
    } else if !r2_records.is_empty() {
        r2_records[0].data().get(&Tag::READ_GROUP).and_then(|v| {
            match v {
                Value::String(s) => Some(s.to_string()),
                _ => None,
            }
        })
    } else {
        None
    };
    
    // Extract extra tag (try first record of R1 as representative)
    let mut extra_tag = None;
    if !r1_records.is_empty() {
        extra_tag = get_tag_value(r1_records[0], &config.extract_tag);
    } else if !r2_records.is_empty() {
        extra_tag = get_tag_value(r2_records[0], &config.extract_tag);
    }
    
    // Check if we have both R1 and R2
    if r1_records.is_empty() || r2_records.is_empty() {
        stats.pair_type_nn += 1;
        // Create fallback using first available records
        let r1_segkeys: Vec<SegKey> = r1_records.iter().map(|r| record_to_segkey(r, config)).collect();
        let r2_segkeys: Vec<SegKey> = r2_records.iter().map(|r| record_to_segkey(r, config)).collect();
        return create_fallback_pair(&read_id, &r1_segkeys, &r2_segkeys, rg, extra_tag);
    }
    
    if r1_records.len() + r2_records.len() > 3 {
        stats.many_hits += 1;
        stats.pair_type_mm += 1;
    }
    
    // Determine pair type based on MapQ
    let r1_0 = record_to_segkey(r1_records[0], config);
    let r2_0 = record_to_segkey(r2_records[0], config);
    let q1_good = r1_0.is_mapped && r1_0.mapq >= config.min_mapq;
    let q2_good = r2_0.is_mapped && r2_0.mapq >= config.min_mapq;
    
    let mut pair_type = if !q1_good && !q2_good {
        "NN"
    } else if !q1_good || !q2_good {
        "NU"
    } else {
        "UU"
    };
    
    // Parse segments and determine category
    let category = if r1_records.len() == 1 && r2_records.len() == 1 {
        0 // 1+1
    } else if r1_records.len() == 1 && r2_records.len() == 2 {
        1 // 1+2
    } else if r1_records.len() == 2 && r2_records.len() == 1 {
        2 // 2+1
    } else {
        // Unsupported configuration
        stats.many_hits += 1;
        return None;
    };
    
    // Parse CIGAR strings from records
    let s1 = parse_segment_from_record(r1_records[0]);
    let s2 = if category == 0 || category == 1 {
        parse_segment_from_record(r2_records[0])
    } else {
        // category == 2, s2 is from second R1 record
        if r1_records.len() > 1 {
            parse_segment_from_record(r1_records[1])
        } else {
            Segment::new()
        }
    };
    let s3 = if category == 1 {
        // s3 is from second R2 record
        if r2_records.len() > 1 {
            parse_segment_from_record(r2_records[1])
        } else {
            Segment::new()
        }
    } else if category == 2 {
        parse_segment_from_record(r2_records[0])
    } else {
        Segment::new()
    };
    
    // Check integrity
    let r1_0 = record_to_segkey(r1_records[0], config);
    let r2_0 = record_to_segkey(r2_records[0], config);
    
    if category == 0 {
        if r1_0.is_mapped && !check_integrity_1_seg(&s1, config) {
            stats.low_map += 1;
            pair_type = "NN";
        }
        if r2_0.is_mapped && !check_integrity_1_seg(&s2, config) {
            stats.low_map += 1;
            pair_type = "NN";
        }
        if r1_0.is_mapped && r2_0.is_mapped && s1.seg_cnt + s2.seg_cnt > 3 {
            stats.many_hits += 1;
            pair_type = "MM";
            stats.pair_type_mm += 1;
        }
    } else if category == 1 {
        if r1_0.is_mapped && !check_integrity_1_seg(&s1, config) {
            stats.low_map += 1;
            pair_type = "NN";
        }
        let r2_1 = if r2_records.len() > 1 {
            record_to_segkey(r2_records[1], config)
        } else {
            record_to_segkey(r2_records[0], config)
        };
        if r2_0.is_mapped && r2_1.is_mapped && !check_integrity_2_seg(&s2, &s3, config) {
            stats.low_map += 1;
            pair_type = "NN";
        }
        if r1_0.is_mapped && r2_0.is_mapped && r2_1.is_mapped
            && (s1.seg_cnt != 1 || s2.seg_cnt != 1 || s3.seg_cnt != 1)
        {
            stats.many_hits += 1;
            pair_type = "MM";
            stats.pair_type_mm += 1;
        }
    } else {
        // category == 2
        let r1_1 = if r1_records.len() > 1 {
            record_to_segkey(r1_records[1], config)
        } else {
            record_to_segkey(r1_records[0], config)
        };
        if r1_0.is_mapped && r1_1.is_mapped && !check_integrity_2_seg(&s1, &s2, config) {
            stats.low_map += 1;
            pair_type = "NN";
        }
        if r2_0.is_mapped && !check_integrity_1_seg(&s3, config) {
            stats.low_map += 1;
            pair_type = "NN";
        }
        if r1_0.is_mapped && r1_1.is_mapped && r2_0.is_mapped
            && (s1.seg_cnt != 1 || s2.seg_cnt != 1 || s3.seg_cnt != 1)
        {
            stats.many_hits += 1;
            pair_type = "MM";
            stats.pair_type_mm += 1;
        }
    }
    
    // Resolve positions based on category
    let result = match category {
        0 => resolve_category_0(&r1_records, &r2_records, &s1, &s2, config, stats, &mut pair_type),
        1 => resolve_category_1(&r1_records, &r2_records, &s1, &s2, &s3, config, stats, &mut pair_type),
        2 => resolve_category_2(&r1_records, &r2_records, &s1, &s2, &s3, config, stats, &mut pair_type),
        _ => return None,
    };
    
    let (chr1, pos1, strand1, start1, end1, chr2, pos2, strand2, start2, end2, q1, q2) = match result {
        Some(r) => r,
        None => return None,
    };
    
    // Calculate pair MapQ
    let pair_mapq = q1.min(q2);
    stats.sum_mapq += pair_mapq as u64;
    stats.valid_pairs += 1;
    
    // Update pair type counters
    match pair_type {
        "UU" => stats.pair_type_uu += 1,
        "NU" => stats.pair_type_nu += 1,
        "NN" => stats.pair_type_nn += 1,
        "MM" => {}, // Already counted
        "UR" => stats.pair_type_ur += 1,
        "SC" => stats.pair_type_sc += 1,
        _ => {}
    }
    
    // Extract RG tag from R1 records (already extracted above)
    // Extract extra tag if not already done
    if extra_tag.is_none() {
        extra_tag = get_tag_value(r1_records[0], &config.extract_tag);
    }
    
    // Determine final order and calculate distance
    let (final_chr1, final_pos1, final_strand1, final_start1, final_end1, final_chr2, final_pos2, final_strand2, final_start2, final_end2, final_mapq1, final_mapq2) =
        if chr1 < chr2 || (chr1 == chr2 && pos1 < pos2) {
            if chr1 == chr2 {
                let dist = if pos2 > pos1 { pos2 - pos1 } else { 0 };
                if dist <= config.max_self_circle_dist {
                    stats.self_circle += 1;
                    return None; // Skip self-circles
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
            (chr1, pos1, strand1, start1, end1, chr2, pos2, strand2, start2, end2, q1, q2)
        } else {
            if chr1 == chr2 {
                let dist = if pos1 > pos2 { pos1 - pos2 } else { 0 };
                if dist <= config.max_self_circle_dist {
                    stats.self_circle += 1;
                    return None; // Skip self-circles
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
            (chr2, pos2, strand2, start2, end2, chr1, pos1, strand1, start1, end1, q2, q1)
        };
    
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
        
        let is_mapped1 = final_chr1 != "!";
        let is_mapped2 = final_chr2 != "!";
        
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

fn record_to_segkey(record: &RecordBuf, config: &Config) -> SegKey {
    let flags = record.flags();
    let mapq = record.mapping_quality().map(|q| q.get()).unwrap_or(0);
    let is_mapped = !flags.is_unmapped()
        && record.reference_sequence_id().is_some()
        && !record.cigar().as_ref().is_empty();
    let chr = if is_mapped {
        config.get_chrom_name(record.reference_sequence_id())
    } else {
        "!".to_string()
    };
    let pos = record.alignment_start().map(|p| p.get()).unwrap_or(0);
    let strand = if flags.is_reverse_complemented() { '-' } else { '+' };
    
    SegKey {
        chr,
        pos,
        strand,
        cigar: String::new(),
        mapq,
        is_mapped,
    }
}

fn parse_segment_from_record(record: &RecordBuf) -> Segment {
    if record.flags().is_unmapped() || record.reference_sequence_id().is_none() {
        return Segment::new();
    }
    
    let start = record.alignment_start().map(|p| p.get()).unwrap_or(0);
    match cigar_to_segment(record.cigar(), start) {
        Ok(seg) => seg,
        Err(_) => Segment::new(),
    }
}

// Helper to extract a specific tag from a record
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

fn create_fallback_pair(
    read_id: &str,
    r1_segkeys: &[SegKey],
    r2_segkeys: &[SegKey],
    rg: Option<String>,
    extra_tag: Option<String>,
) -> Option<PairOutput> {
    let (chr1, pos1, strand1, q1) = if !r1_segkeys.is_empty() {
        (
            r1_segkeys[0].chr.clone(),
            r1_segkeys[0].pos,
            r1_segkeys[0].strand,
            r1_segkeys[0].mapq,
        )
    } else {
        ("!".to_string(), 0, '-', 0)
    };
    
    let (chr2, pos2, strand2, q2) = if !r2_segkeys.is_empty() {
        (
            r2_segkeys[0].chr.clone(),
            r2_segkeys[0].pos,
            r2_segkeys[0].strand,
            r2_segkeys[0].mapq,
        )
    } else {
        ("!".to_string(), 0, '-', 0)
    };
    
    Some(PairOutput {
        read_id: read_id.to_string(),
        chr1,
        pos1,
        chr2,
        pos2,
        strand1,
        strand2,
        pair_type: "NN".to_string(),
        mapq1: q1,
        mapq2: q2,
        rg,
        start1: 0,
        end1: 0,
        start2: 0,
        end2: 0,
        extra_tag,
        hic_type: None, // Fallback pairs don't have fragment information
    })
}

// Updated return type: (chr1, pos1, strand1, start1, end1, chr2, pos2, strand2, start2, end2, q1, q2)
type ResolveResult = Option<(String, usize, char, usize, usize, String, usize, char, usize, usize, u8, u8)>;

fn resolve_category_0(
    r1_records: &[&RecordBuf],
    r2_records: &[&RecordBuf],
    s1: &Segment,
    s2: &Segment,
    config: &Config,
    _stats: &mut Statistics,
    _pair_type: &mut &str,
) -> ResolveResult {
    let r1_0 = record_to_segkey(r1_records[0], config);
    let r2_0 = record_to_segkey(r2_records[0], config);
    let strand1 = r1_0.strand;
    let strand2 = r2_0.strand;
    let chr1 = r1_0.chr.clone();
    let chr2 = r2_0.chr.clone();
    
    let pos1 = if s1.seg_cnt == 1 && s1.left.len() > 0 {
        if strand1 == '+' {
            s1.left[0]
        } else {
            s1.right[0]
        }
    } else {
        r1_0.pos
    };
    
    let pos2 = if s2.seg_cnt == 1 && s2.left.len() > 0 {
        if strand2 == '+' {
            s2.left[0]
        } else {
            s2.right[0]
        }
    } else {
        r2_0.pos
    };
    
    let (start1, end1) = s1.span();
    let (start2, end2) = s2.span();
    
    Some((chr1, pos1, strand1, start1, end1, chr2, pos2, strand2, start2, end2, r1_0.mapq, r2_0.mapq))
}

fn resolve_category_1(
    r1_records: &[&RecordBuf],
    r2_records: &[&RecordBuf],
    s1: &Segment,
    s2: &Segment,
    s3: &Segment,
    config: &Config,
    stats: &mut Statistics,
    pair_type: &mut &str,
) -> ResolveResult {
    let r1_0 = record_to_segkey(r1_records[0], config);
    let r2_0 = record_to_segkey(r2_records[0], config);
    let r2_1 = if r2_records.len() > 1 {
        record_to_segkey(r2_records[1], config)
    } else {
        record_to_segkey(r2_records[0], config)
    };
    let strand1 = r1_0.strand;
    let chr1 = r1_0.chr.clone();
    
    let pos1 = if s1.left.len() > 0 {
        if strand1 == '+' {
            s1.left[0]
        } else {
            s1.right[0]
        }
    } else {
        r1_0.pos
    };
    
    let (start1, end1) = s1.span();
    
    // Try to pair s1 with s2 or s3
    let mate = determine_mate_1_2(&r1_0, &r2_0, &r2_1, s1, s2, s3, config);
    
    if mate == 0 {
        stats.unpaired += 1;
        *pair_type = "UR";
        stats.pair_type_ur += 1;
        
        // Fallback to primary alignments
        let chr2 = r2_0.chr.clone();
        let strand2 = r2_0.strand;
        let pos2 = if s2.left.len() > 0 {
            if strand2 == '+' {
                s2.left[0]
            } else {
                s2.right[0]
            }
        } else {
            r2_0.pos
        };
        let (start2, end2) = s2.span();
        
        return Some((chr1, pos1, strand1, start1, end1, chr2, pos2, strand2, start2, end2, r1_0.mapq, r2_0.mapq));
    }
    
    let (chr2, pos2, strand2, start2, end2, q2) = if mate == 2 {
        // s1 paired with s2, extract from s3
        let chr = r2_1.chr.clone();
        let strand = r2_1.strand;
        let pos = if s3.left_clip > s3.right_clip && s3.right.len() > 0 {
            s3.right[0]
        } else if s3.left.len() > 0 {
            s3.left[0]
        } else {
            r2_1.pos
        };
        let (s, e) = s3.span();
        (chr, pos, strand, s, e, r2_0.mapq)
    } else {
        // mate == 3: s1 paired with s3, extract from s2
        let chr = r2_0.chr.clone();
        let strand = r2_0.strand;
        let pos = if s2.left_clip > s2.right_clip && s2.right.len() > 0 {
            s2.right[0]
        } else if s2.left.len() > 0 {
            s2.left[0]
        } else {
            r2_0.pos
        };
        let (s, e) = s2.span();
        (chr, pos, strand, s, e, r2_1.mapq)
    };
    
    Some((chr1, pos1, strand1, start1, end1, chr2, pos2, strand2, start2, end2, r1_0.mapq, q2))
}

fn resolve_category_2(
    r1_records: &[&RecordBuf],
    r2_records: &[&RecordBuf],
    s1: &Segment,
    s2: &Segment,
    s3: &Segment,
    config: &Config,
    stats: &mut Statistics,
    pair_type: &mut &str,
) -> ResolveResult {
    let r1_0 = record_to_segkey(r1_records[0], config);
    let r1_1 = if r1_records.len() > 1 {
        record_to_segkey(r1_records[1], config)
    } else {
        record_to_segkey(r1_records[0], config)
    };
    let r2_0 = record_to_segkey(r2_records[0], config);
    let strand2 = r2_0.strand;
    let chr2 = r2_0.chr.clone();
    
    let pos2 = if s3.left.len() > 0 {
        if strand2 == '+' {
            s3.left[0]
        } else {
            s3.right[0]
        }
    } else {
        r2_0.pos
    };
    
    let (start2, end2) = s3.span();
    
    // Try to pair s3 with s1 or s2
    let mate = determine_mate_2_1(&r1_0, &r1_1, &r2_0, s1, s2, s3, config);
    
    if mate == 0 {
        stats.unpaired += 1;
        *pair_type = "UR";
        stats.pair_type_ur += 1;
        
        // Fallback to primary alignments
        let chr1 = r1_0.chr.clone();
        let strand1 = r1_0.strand;
        let pos1 = if s1.left.len() > 0 {
            if strand1 == '+' {
                s1.left[0]
            } else {
                s1.right[0]
            }
        } else {
            r1_0.pos
        };
        let (start1, end1) = s1.span();
        
        return Some((chr1, pos1, strand1, start1, end1, chr2, pos2, strand2, start2, end2, r1_0.mapq, r2_0.mapq));
    }
    
    let (chr1, pos1, strand1, start1, end1, q1) = if mate == 1 {
        // s3 paired with s1, extract from s2
        let chr = r1_1.chr.clone();
        let strand = r1_1.strand;
        let pos = if s2.left_clip > s2.right_clip && s2.right.len() > 0 {
            s2.right[0]
        } else if s2.left.len() > 0 {
            s2.left[0]
        } else {
            r1_1.pos
        };
        let (s, e) = s2.span();
        (chr, pos, strand, s, e, r1_0.mapq)
    } else {
        // mate == 2: s3 paired with s2, extract from s1
        let chr = r1_0.chr.clone();
        let strand = r1_0.strand;
        let pos = if s1.left_clip > s1.right_clip && s1.right.len() > 0 {
            s1.right[0]
        } else if s1.left.len() > 0 {
            s1.left[0]
        } else {
            r1_0.pos
        };
        let (s, e) = s1.span();
        (chr, pos, strand, s, e, r1_1.mapq)
    };
    
    Some((chr1, pos1, strand1, start1, end1, chr2, pos2, strand2, start2, end2, q1, r2_0.mapq))
}

fn determine_mate_1_2(
    r1_0: &SegKey,
    r2_0: &SegKey,
    r2_1: &SegKey,
    s1: &Segment,
    s2: &Segment,
    s3: &Segment,
    config: &Config,
) -> usize {
    let strand1 = r1_0.strand;
    
    // Check s1 vs s2
    if strand1 == '+' {
        if r2_0.strand == '-'
            && r1_0.chr == r2_0.chr
            && s1.left.len() > 0
            && s2.left.len() > 0
            && s1.left[0] < s2.left[0]
            && s2.right[0].saturating_sub(s1.left[0]) <= config.max_pair_dist
        {
            return 2;
        }
    } else {
        if r2_0.strand == '+'
            && r1_0.chr == r2_0.chr
            && s2.left.len() > 0
            && s1.left.len() > 0
            && s2.left[0] < s1.left[0]
            && s1.right[0].saturating_sub(s2.left[0]) <= config.max_pair_dist
        {
            return 2;
        }
    }
    
    // Check s1 vs s3
    if strand1 == '+' {
        if r2_1.strand == '-'
            && r1_0.chr == r2_1.chr
            && s1.left.len() > 0
            && s3.left.len() > 0
            && s1.left[0] < s3.left[0]
            && s3.right[0].saturating_sub(s1.left[0]) <= config.max_pair_dist
        {
            return 3;
        }
    } else {
        if r2_1.strand == '+'
            && r1_0.chr == r2_1.chr
            && s3.left.len() > 0
            && s1.left.len() > 0
            && s3.left[0] < s1.left[0]
            && s1.right[0].saturating_sub(s3.left[0]) <= config.max_pair_dist
        {
            return 3;
        }
    }
    
    0
}

fn determine_mate_2_1(
    r1_0: &SegKey,
    r1_1: &SegKey,
    r2_0: &SegKey,
    s1: &Segment,
    s2: &Segment,
    s3: &Segment,
    config: &Config,
) -> usize {
    let strand2 = r2_0.strand;
    
    // Check s3 vs s1
    if strand2 == '+' {
        if r1_0.strand == '-'
            && r2_0.chr == r1_0.chr
            && s3.left.len() > 0
            && s1.left.len() > 0
            && s3.left[0] < s1.left[0]
            && s1.right[0].saturating_sub(s3.left[0]) <= config.max_pair_dist
        {
            return 1;
        }
    } else {
        if r1_0.strand == '+'
            && r2_0.chr == r1_0.chr
            && s1.left.len() > 0
            && s3.left.len() > 0
            && s1.left[0] < s3.left[0]
            && s3.right[0].saturating_sub(s1.left[0]) <= config.max_pair_dist
        {
            return 1;
        }
    }
    
    // Check s3 vs s2
    if strand2 == '+' {
        if r1_1.strand == '-'
            && r2_0.chr == r1_1.chr
            && s3.left.len() > 0
            && s2.left.len() > 0
            && s3.left[0] < s2.left[0]
            && s2.right[0].saturating_sub(s3.left[0]) <= config.max_pair_dist
        {
            return 2;
        }
    } else {
        if r1_1.strand == '+'
            && r2_0.chr == r1_1.chr
            && s2.left.len() > 0
            && s3.left.len() > 0
            && s2.left[0] < s3.left[0]
            && s3.right[0].saturating_sub(s2.left[0]) <= config.max_pair_dist
        {
            return 2;
        }
    }
    
    0
}
