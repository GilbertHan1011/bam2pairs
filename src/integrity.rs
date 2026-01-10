use crate::config::Config;
use crate::segment::Segment;

/// Check integrity of a single segment
/// Port of check_integrity_1_seg from pairutil.h (lines 190-198)
pub fn check_integrity_1_seg(segment: &Segment, config: &Config) -> bool {
    let mut total = segment.mappable;
    
    // Add clips if they exceed minimum size
    if segment.left_clip > config.min_clip_size {
        total += segment.left_clip;
    }
    if segment.right_clip > config.min_clip_size {
        total += segment.right_clip;
    }
    
    // Check if mappable bases meet minimum ratio
    segment.mappable >= (total as f32 * config.min_mapped_ratio) as usize
}

/// Check integrity of two segments
/// Port of check_integrity_2_seg from pairutil.h (lines 200-218)
pub fn check_integrity_2_seg(seg1: &Segment, seg2: &Segment, config: &Config) -> bool {
    let mut total_1 = seg1.mappable;
    if seg1.left_clip > config.min_clip_size {
        total_1 += seg1.left_clip;
    }
    if seg1.right_clip > config.min_clip_size {
        total_1 += seg1.right_clip;
    }
    
    let mut total_2 = seg2.mappable;
    if seg2.left_clip > config.min_clip_size {
        total_2 += seg2.left_clip;
    }
    if seg2.right_clip > config.min_clip_size {
        total_2 += seg2.right_clip;
    }
    
    let total_mappable = seg1.mappable + seg2.mappable;
    let max_total = total_1.max(total_2);
    
    total_mappable >= (max_total as f32 * config.min_mapped_ratio) as usize
}

#[cfg(test)]
mod tests {
    use super::*;
    
    fn create_segment(mappable: usize, left_clip: usize, right_clip: usize) -> Segment {
        let mut seg = Segment::new();
        seg.mappable = mappable;
        seg.left_clip = left_clip;
        seg.right_clip = right_clip;
        seg.seg_cnt = 1;
        seg.left.push(1000);
        seg.right.push(1000 + mappable - 1);
        seg
    }
    
    #[test]
    fn test_integrity_1_seg_pass() {
        let config = Config::default(); // min_mapped_ratio = 0.5
        let seg = create_segment(100, 5, 5); // 100 mappable, 5+5 clips (ignored as < 20)
        assert!(check_integrity_1_seg(&seg, &config));
    }
    
    #[test]
    fn test_integrity_1_seg_fail() {
        let config = Config::default();
        let seg = create_segment(40, 30, 30); // 40 mappable, 60 total with clips
        assert!(!check_integrity_1_seg(&seg, &config)); // 40 < 60 * 0.5
    }
    
    #[test]
    fn test_integrity_1_seg_with_large_clips() {
        let config = Config::default();
        let seg = create_segment(60, 25, 25); // 60 mappable, 110 total (60+25+25)
        assert!(check_integrity_1_seg(&seg, &config)); // 60 >= 110 * 0.5 = 55
    }
    
    #[test]
    fn test_integrity_2_seg_pass() {
        let config = Config::default();
        let seg1 = create_segment(50, 5, 5);
        let seg2 = create_segment(50, 5, 5);
        assert!(check_integrity_2_seg(&seg1, &seg2, &config)); // 100 >= 100 * 0.5
    }
    
    #[test]
    fn test_integrity_2_seg_fail() {
        let config = Config::default();
        let seg1 = create_segment(30, 25, 25); // total_1 = 80 (30 + 25 + 25), mappable = 30
        let seg2 = create_segment(20, 5, 5);   // total_2 = 30 (20 + 5 + 5, but 5 < 20 so ignored), mappable = 20
        // total_mappable = 50, max_total = 80, required = 80 * 0.5 = 40, so 50 >= 40 should pass
        // Let's make it fail: seg1 = 20 mappable, seg2 = 10 mappable, total = 30 < 40
        let seg1_fail = create_segment(20, 25, 25); // total_1 = 70, mappable = 20
        let seg2_fail = create_segment(10, 5, 5);   // total_2 = 20, mappable = 10
        assert!(!check_integrity_2_seg(&seg1_fail, &seg2_fail, &config)); // 30 < 70 * 0.5 = 35
    }
}
