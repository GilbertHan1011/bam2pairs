use noodles::sam::alignment::RecordBuf;
use std::io::BufRead;

/// Group BAM records by read name
/// Returns an iterator of Vec<RecordBuf> where each Vec contains all records with the same name
pub struct ReadGrouper<R: BufRead> {
    reader: noodles::bam::io::Reader<R>,
    header: noodles::sam::Header,
    current_group: Vec<RecordBuf>,
    last_name: Option<Vec<u8>>,
    finished: bool,
}

impl<R: BufRead> ReadGrouper<R> {
    pub fn new(reader: noodles::bam::io::Reader<R>, header: noodles::sam::Header) -> Self {
        Self {
            reader,
            header,
            current_group: Vec::new(),
            last_name: None,
            finished: false,
        }
    }
    
    /// Get the next group of records with the same read name
    pub fn next_group(&mut self) -> Option<Vec<RecordBuf>> {
        if self.finished {
            return None;
        }
        
        loop {
            let mut record = RecordBuf::default();
            match self.reader.read_record_buf(&self.header, &mut record) {
                Ok(0) => {
                    // EOF reached
                    self.finished = true;
                    if !self.current_group.is_empty() {
                        let group = std::mem::take(&mut self.current_group);
                        return Some(group);
                    }
                    return None;
                }
                Ok(_) => {
                    let current_name = record.name().map(|n| n.to_vec());
                    
                    // Skip secondary alignments
                    if record.flags().is_secondary() {
                        continue;
                    }
                    
                    match (&self.last_name, &current_name) {
                        (Some(last), Some(curr)) if last == curr => {
                            // Same read name, add to current group
                            self.current_group.push(record);
                        }
                        (Some(_), Some(curr)) => {
                            // Different read name, return current group and start new one
                            let group = std::mem::take(&mut self.current_group);
                            self.current_group.push(record);
                            self.last_name = Some(curr.clone());
                            return Some(group);
                        }
                        (None, Some(curr)) => {
                            // First record
                            self.current_group.push(record);
                            self.last_name = Some(curr.clone());
                        }
                        _ => {
                            // Unnamed record, skip
                            continue;
                        }
                    }
                }
                Err(_) => {
                    // Error reading, finish
                    self.finished = true;
                    if !self.current_group.is_empty() {
                        let group = std::mem::take(&mut self.current_group);
                        return Some(group);
                    }
                    return None;
                }
            }
        }
    }
}

/// Collect all groups into a Vec for parallel processing
pub fn collect_groups<R: BufRead>(reader: noodles::bam::io::Reader<R>, header: noodles::sam::Header) -> Vec<Vec<RecordBuf>> {
    let mut grouper = ReadGrouper::new(reader, header);
    let mut groups = Vec::new();
    
    while let Some(group) = grouper.next_group() {
        groups.push(group);
    }
    
    groups
}
