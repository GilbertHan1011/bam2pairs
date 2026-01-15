use crate::cigar::cigar_to_segment;
use crate::config::Config;
use crate::flash::PairOutput;
use crate::hic_classifier::{classify_interaction, ReadInfo as HicReadInfo};
use crate::integrity::{check_integrity_1_seg, check_integrity_2_seg};
use crate::segment::{Segment, SegKey, Statistics};
use crate::fragment::Fragment;
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
    
    // Resolve positions based on category and precompute fragments if available
    let result = match category {
        0 => resolve_category_0(&r1_records, &r2_records, &s1, &s2, config, stats, &mut pair_type),
        1 => resolve_category_1(&r1_records, &r2_records, &s1, &s2, &s3, config, stats, &mut pair_type),
        2 => resolve_category_2(&r1_records, &r2_records, &s1, &s2, &s3, config, stats, &mut pair_type),
        _ => return None,
    };
    
    let (chr1, pos1, strand1, start1, end1, chr2, pos2, strand2, start2, end2, q1, q2, frag1_opt, frag2_opt) =
        match result {
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
    let (final_chr1, final_pos1, final_strand1, final_start1, final_end1, final_chr2, final_pos2, final_strand2, final_start2, final_end2, final_mapq1, final_mapq2, final_frag1, final_frag2) =
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
            (chr1, pos1, strand1, start1, end1, chr2, pos2, strand2, start2, end2, q1, q2, frag1_opt.clone(), frag2_opt.clone())
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
            (chr2, pos2, strand2, start2, end2, chr1, pos1, strand1, start1, end1, q2, q1, frag2_opt.clone(), frag1_opt.clone())
        };
    
    // Classify HiC interaction if fragment map is available
    let hic_type = if let Some(ref frag_map) = config.fragment_map {
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
            final_frag1.as_ref(),
            final_frag2.as_ref(),
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

// Updated return type:
// (chr1, pos1, strand1, start1, end1, chr2, pos2, strand2, start2, end2, q1, q2, frag1, frag2)
type ResolveResult = Option<(
    String,
    usize,
    char,
    usize,
    usize,
    String,
    usize,
    char,
    usize,
    usize,
    u8,
    u8,
    Option<Fragment>,
    Option<Fragment>,
)>;

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

    // Precompute fragments if fragment map is available
    let (frag1, frag2) = if config.fragment_map.is_some() {
        let f1 = get_fragment(config, &r1_0, s1);
        let f2 = get_fragment(config, &r2_0, s2);
        (f1, f2)
    } else {
        (None, None)
    };
    
    Some((
        chr1,
        pos1,
        strand1,
        start1,
        end1,
        chr2,
        pos2,
        strand2,
        start2,
        end2,
        r1_0.mapq,
        r2_0.mapq,
        frag1,
        frag2,
    ))
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

    // Precompute fragments if fragment map is available
    let (f1, f2, f3) = if config.fragment_map.is_some() {
        (
            get_fragment(config, &r1_0, s1),
            get_fragment(config, &r2_0, s2),
            get_fragment(config, &r2_1, s3),
        )
    } else {
        (None, None, None)
    };
    
    // Try to pair s1 with s2 or s3 using updated logic
    let mate = determine_mate_1_2(&r1_0, &r2_0, &r2_1, s1, s2, s3, &f1, &f2, &f3, config);
    
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

        // Use fragments corresponding to s1 and s2
        let frag1 = f1.clone();
        let frag2 = f2.clone();
        
        return Some((
            chr1,
            pos1,
            strand1,
            start1,
            end1,
            chr2,
            pos2,
            strand2,
            start2,
            end2,
            r1_0.mapq,
            r2_0.mapq,
            frag1,
            frag2,
        ));
    }
    
    let (chr2, pos2, strand2, start2, end2, q2, frag1, frag2) = if mate == 2 {
        // mate == 2 means s1 and s2 are local; output contact s1-s3
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

        // final fragments correspond to s1 and s3
        (chr, pos, strand, s, e, r2_1.mapq, f1.clone(), f3.clone())
    } else {
        // mate == 3 means s1 and s3 are local; output contact s1-s2
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

        // final fragments correspond to s1 and s2
        (chr, pos, strand, s, e, r2_0.mapq, f1.clone(), f2.clone())
    };
    
    Some((
        chr1,
        pos1,
        strand1,
        start1,
        end1,
        chr2,
        pos2,
        strand2,
        start2,
        end2,
        r1_0.mapq,
        q2,
        frag1,
        frag2,
    ))
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

    // Precompute fragments if fragment map is available
    let (f1, f2, f3) = if config.fragment_map.is_some() {
        (
            get_fragment(config, &r1_0, s1),
            get_fragment(config, &r1_1, s2),
            get_fragment(config, &r2_0, s3),
        )
    } else {
        (None, None, None)
    };
    
    // Try to pair s3 with s1 or s2 using updated logic
    let mate = determine_mate_2_1(&r1_0, &r1_1, &r2_0, s1, s2, s3, &f1, &f2, &f3, config);
    
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

        // Use fragments corresponding to s3 and s1
        let frag1 = f1.clone();
        let frag2 = f3.clone();
        
        return Some((
            chr1,
            pos1,
            strand1,
            start1,
            end1,
            chr2,
            pos2,
            strand2,
            start2,
            end2,
            r1_0.mapq,
            r2_0.mapq,
            frag1,
            frag2,
        ));
    }
    
    let (chr1, pos1, strand1, start1, end1, q1, frag1, frag2) = if mate == 1 {
        // mate == 1 means s3 and s1 are local; output contact s3-s2
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

        // final fragments correspond to s2 and s3
        (chr, pos, strand, s, e, r1_1.mapq, f2.clone(), f3.clone())
    } else {
        // mate == 2 means s3 and s2 are local; output contact s3-s1
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

        // final fragments correspond to s1 and s3
        (chr, pos, strand, s, e, r1_0.mapq, f1.clone(), f3.clone())
    };
    
    Some((
        chr1,
        pos1,
        strand1,
        start1,
        end1,
        chr2,
        pos2,
        strand2,
        start2,
        end2,
        q1,
        r2_0.mapq,
        frag1,
        frag2,
    ))
}

fn determine_mate_1_2(
    r1_0: &SegKey,
    r2_0: &SegKey,
    r2_1: &SegKey,
    s1: &Segment,
    s2: &Segment,
    s3: &Segment,
    f1: &Option<Fragment>,
    f2: &Option<Fragment>,
    f3: &Option<Fragment>,
    config: &Config,
) -> usize {
    // 1. Restriction Fragment Check (Biological Truth)
    if config.fragment_map.is_some() {
        // If s1 and s2 are same/adjacent fragment -> s2 is the local mate
        if are_same_fragment(f1, f2) || are_contiguous(f1, f2) {
            return 2;
        }
        // If s1 and s3 are same/adjacent fragment -> s3 is the local mate
        if are_same_fragment(f1, f3) || are_contiguous(f1, f3) {
            return 3;
        }
    }

    // 2. Collinearity Check (Physical Geometry Fallback)
    let strand1 = r1_0.strand;
    
    // Check s1 vs s2 (Are they inward facing and close?)
    let s1_s2_local = if strand1 == '+' {
        r2_0.strand == '-'
            && r1_0.chr == r2_0.chr
            && s1.left.first().unwrap_or(&0) < s2.left.first().unwrap_or(&0)
            && s2
                .right
                .first()
                .unwrap_or(&0)
                .saturating_sub(*s1.left.first().unwrap_or(&0))
                <= config.max_pair_dist
    } else {
        r2_0.strand == '+'
            && r1_0.chr == r2_0.chr
            && s2.left.first().unwrap_or(&0) < s1.left.first().unwrap_or(&0)
            && s1
                .right
                .first()
                .unwrap_or(&0)
                .saturating_sub(*s2.left.first().unwrap_or(&0))
                <= config.max_pair_dist
    };

    if s1_s2_local {
        return 2;
    }

    // Check s1 vs s3
    let s1_s3_local = if strand1 == '+' {
        r2_1.strand == '-'
            && r1_0.chr == r2_1.chr
            && s1.left.first().unwrap_or(&0) < s3.left.first().unwrap_or(&0)
            && s3
                .right
                .first()
                .unwrap_or(&0)
                .saturating_sub(*s1.left.first().unwrap_or(&0))
                <= config.max_pair_dist
    } else {
        r2_1.strand == '+'
            && r1_0.chr == r2_1.chr
            && s3.left.first().unwrap_or(&0) < s1.left.first().unwrap_or(&0)
            && s1
                .right
                .first()
                .unwrap_or(&0)
                .saturating_sub(*s3.left.first().unwrap_or(&0))
                <= config.max_pair_dist
    };

    if s1_s3_local {
        return 3;
    }

    // 3. MapQ Check (Statistical Ambiguity Resolution)
    // If neither looks like a local pair (Trans-Trans or Long-Range),
    // we assume the High MapQ alignment is the valid CONTACT,
    // and the Low MapQ alignment is the noise/artifact.
    //
    // Note: 'mate' is the index of the segment we *don't* want in the output pair.
    
    if r2_0.mapq >= r2_1.mapq {
        // s2 is better quality. We want output to be s1-s2.
        // So we designate s3 as the "mate" (to be hidden/extracted away).
        3
    } else {
        // s3 is better quality. We want output to be s1-s3.
        // So we designate s2 as the "mate".
        2
    }
}

fn determine_mate_2_1(
    r1_0: &SegKey,
    r1_1: &SegKey,
    r2_0: &SegKey,
    s1: &Segment,
    s2: &Segment,
    s3: &Segment,
    f1: &Option<Fragment>,
    f2: &Option<Fragment>,
    f3: &Option<Fragment>,
    config: &Config,
) -> usize {
    // 1. Restriction Fragment Check
    if config.fragment_map.is_some() {
        // s3 vs s1
        if are_same_fragment(f3, f1) || are_contiguous(f3, f1) {
            return 1;
        }
        // s3 vs s2
        if are_same_fragment(f3, f2) || are_contiguous(f3, f2) {
            return 2;
        }
    }

    // 2. Collinearity Check
    let strand2 = r2_0.strand;
    
    // Check s3 vs s1
    let s3_s1_local = if strand2 == '+' {
        r1_0.strand == '-'
            && r2_0.chr == r1_0.chr
            && s3.left.first().unwrap_or(&0) < s1.left.first().unwrap_or(&0)
            && s1
                .right
                .first()
                .unwrap_or(&0)
                .saturating_sub(*s3.left.first().unwrap_or(&0))
                <= config.max_pair_dist
    } else {
        r1_0.strand == '+'
            && r2_0.chr == r1_0.chr
            && s1.left.first().unwrap_or(&0) < s3.left.first().unwrap_or(&0)
            && s3
                .right
                .first()
                .unwrap_or(&0)
                .saturating_sub(*s1.left.first().unwrap_or(&0))
                <= config.max_pair_dist
    };

    if s3_s1_local {
        return 1;
    }

    // Check s3 vs s2
    let s3_s2_local = if strand2 == '+' {
        r1_1.strand == '-'
            && r2_0.chr == r1_1.chr
            && s3.left.first().unwrap_or(&0) < s2.left.first().unwrap_or(&0)
            && s2
                .right
                .first()
                .unwrap_or(&0)
                .saturating_sub(*s3.left.first().unwrap_or(&0))
                <= config.max_pair_dist
    } else {
        r1_1.strand == '+'
            && r2_0.chr == r1_1.chr
            && s2.left.first().unwrap_or(&0) < s3.left.first().unwrap_or(&0)
            && s3
                .right
                .first()
                .unwrap_or(&0)
                .saturating_sub(*s2.left.first().unwrap_or(&0))
                <= config.max_pair_dist
    };

    if s3_s2_local {
        return 2;
    }

    // 3. MapQ Check
    // Pick the best quality alignment as the valid contact.
    if r1_0.mapq >= r1_1.mapq {
        // s1 is better. Output s3-s1. Hide s2.
        2
    } else {
        // s2 is better. Output s3-s2. Hide s1.
        1
    }
}

// -----------------------------------------------------------
//  NEW HELPERS FOR FRAGMENT TOPOLOGY
// -----------------------------------------------------------

fn get_fragment(config: &Config, key: &SegKey, seg: &Segment) -> Option<Fragment> {
    let frag_map = config.fragment_map.as_ref()?;
    
    // Calculate middle pos similar to classification
    let (start, end) = seg.span();
    let mid = if start > 0 && end > 0 {
        (start + end) / 2
    } else {
        key.pos
    };
    
    if key.chr == "!" {
        return None;
    }
    frag_map.query_fragment(&key.chr, mid)
}

fn are_same_fragment(f1: &Option<Fragment>, f2: &Option<Fragment>) -> bool {
    match (f1, f2) {
        (Some(a), Some(b)) => a.chr == b.chr && a.start == b.start && a.end == b.end,
        _ => false,
    }
}

fn are_contiguous(f1: &Option<Fragment>, f2: &Option<Fragment>) -> bool {
    match (f1, f2) {
        (Some(a), Some(b)) => {
            if a.chr != b.chr {
                return false;
            }
            a.end == b.start || a.start == b.end
        }
        _ => false,
    }
}
