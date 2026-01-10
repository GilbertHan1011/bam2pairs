use crate::segment::Segment;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record_buf::Cigar;

/// Parse CIGAR string into Segment structure
/// Port of cigar2segment from pairutil.h (lines 74-137)
pub fn cigar_to_segment(cigar: &Cigar, start: usize) -> Result<Segment, String> {
    let mut segment = Segment::new();
    
    // Initialize first segment
    let mut index = 0;
    segment.left.push(start);
    segment.right.push(0);
    
    let mut curr_pos = start;
    let ops: Vec<_> = cigar.as_ref().iter().collect();
    
    for (i, op) in ops.iter().enumerate() {
        let kind = op.kind();
        let length = op.len();
        
        match kind {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                // M, = or X: match/mismatch, consume reference and query
                segment.mappable += length;
                curr_pos += length;
                segment.right[index] = curr_pos - 1;
            }
            Kind::SoftClip | Kind::HardClip => {
                // S or H: clipping
                if i == ops.len() - 1 {
                    // Last operation: right clip
                    segment.right_clip = length;
                } else if index == 0 {
                    // First segment: left clip
                    segment.left_clip = length;
                } else {
                    return Err("ERROR: Clip in middle of CIGAR".to_string());
                }
            }
            Kind::Deletion => {
                // D: deletion, move reference cursor
                curr_pos += length;
                segment.right[index] = curr_pos - 1;
            }
            Kind::Insertion => {
                // I: insertion, skip it (doesn't consume reference)
                // Do nothing
            }
            Kind::Skip => {
                // N: intron/skip, initiate new segment
                curr_pos += length;
                index += 1;
                segment.left.push(curr_pos);
                segment.right.push(0);
            }
            Kind::Pad => {
                // P: padding, skip
                // Do nothing
            }
        }
    }
    
    // Double-check
    if segment.right[index] == 0 {
        return Err(format!("ERROR: Something is wrong with CIGAR parsing"));
    }
    
    segment.seg_cnt = index + 1;
    Ok(segment)
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles::sam::alignment::record::cigar::op::Kind;
    use noodles::sam::alignment::record::cigar::Op;
    
    #[test]
    fn test_simple_match() {
        // 100M
        let cigar = Cigar::try_from(vec![Op::new(Kind::Match, 100)]).unwrap();
        let segment = cigar_to_segment(&cigar, 1000).unwrap();
        
        assert_eq!(segment.seg_cnt, 1);
        assert_eq!(segment.mappable, 100);
        assert_eq!(segment.left[0], 1000);
        assert_eq!(segment.right[0], 1099);
        assert_eq!(segment.left_clip, 0);
        assert_eq!(segment.right_clip, 0);
    }
    
    #[test]
    fn test_with_clips() {
        // 5S95M10S
        let cigar = Cigar::try_from(vec![
            Op::new(Kind::SoftClip, 5),
            Op::new(Kind::Match, 95),
            Op::new(Kind::SoftClip, 10),
        ]).unwrap();
        let segment = cigar_to_segment(&cigar, 1000).unwrap();
        
        assert_eq!(segment.seg_cnt, 1);
        assert_eq!(segment.mappable, 95);
        assert_eq!(segment.left_clip, 5);
        assert_eq!(segment.right_clip, 10);
        assert_eq!(segment.left[0], 1000);
        assert_eq!(segment.right[0], 1094);
    }
    
    #[test]
    fn test_with_intron() {
        // 50M200N50M
        let cigar = Cigar::try_from(vec![
            Op::new(Kind::Match, 50),
            Op::new(Kind::Skip, 200),
            Op::new(Kind::Match, 50),
        ]).unwrap();
        let segment = cigar_to_segment(&cigar, 1000).unwrap();
        
        assert_eq!(segment.seg_cnt, 2);
        assert_eq!(segment.mappable, 100);
        assert_eq!(segment.left[0], 1000);
        assert_eq!(segment.right[0], 1049);
        assert_eq!(segment.left[1], 1250);
        assert_eq!(segment.right[1], 1299);
    }
    
    #[test]
    fn test_with_deletion() {
        // 50M10D50M
        let cigar = Cigar::try_from(vec![
            Op::new(Kind::Match, 50),
            Op::new(Kind::Deletion, 10),
            Op::new(Kind::Match, 50),
        ]).unwrap();
        let segment = cigar_to_segment(&cigar, 1000).unwrap();
        
        assert_eq!(segment.seg_cnt, 1);
        assert_eq!(segment.mappable, 100);
        assert_eq!(segment.left[0], 1000);
        assert_eq!(segment.right[0], 1109); // 1000 + 50 + 10 + 50 - 1
    }
    
    #[test]
    fn test_with_insertion() {
        // 50M10I50M
        let cigar = Cigar::try_from(vec![
            Op::new(Kind::Match, 50),
            Op::new(Kind::Insertion, 10),
            Op::new(Kind::Match, 50),
        ]).unwrap();
        let segment = cigar_to_segment(&cigar, 1000).unwrap();
        
        assert_eq!(segment.seg_cnt, 1);
        assert_eq!(segment.mappable, 100);
        assert_eq!(segment.left[0], 1000);
        assert_eq!(segment.right[0], 1099); // Insertion doesn't consume reference
    }
}
