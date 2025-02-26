use std::cmp::min;
use lib_wfa2::affine_wavefront::{AffineWavefronts, AlignmentStatus, HeuristicStrategy};

/// Convert a CIGAR string into tracepoints.
/// 
/// This function segments a CIGAR string into tracepoints where each segment contains
/// at most `max_diff` differences (mismatches or indels). Match operations ('=' and 'M') 
/// don't count as differences. Key features:
/// 
/// - Mismatch operations ('X') can be split across segments if needed
/// - Indels ('I', 'D') are kept intact within a single segment when possible
/// - Long indels exceeding max_diff become their own segments
/// 
/// @param cigar: The CIGAR string to process
/// @param max_diff: Maximum number of differences allowed in each segment
/// @return Vector of tracepoints, each containing (a_segment_length, b_segment_length)
pub fn cigar_to_tracepoints(
    cigar: &str,
    max_diff: usize,
) -> Vec<(usize, usize)> {
    let ops = cigar_str_to_cigar_ops(cigar);
    let mut tracepoints = Vec::new();
    // current segment counters:
    let mut cur_a_len = 0;
    let mut cur_b_len = 0;
    let mut cur_diff = 0;

    for (mut len, op) in ops {
        match op {
            'X' => {
                // X is splittable
                while len > 0 {
                    let remaining = max_diff.saturating_sub(cur_diff);
                    let step = min(len, remaining);
                    cur_a_len += step;
                    cur_b_len += step;
                    cur_diff += step;
                    len -= step;
                    if cur_diff == max_diff {
                        tracepoints.push((cur_a_len, cur_b_len));
                        cur_a_len = 0;
                        cur_b_len = 0;
                        cur_diff = 0;
                    }
                }
            },
            'I' | 'D' => {
                // For indels, which are unsplittable, try to incorporate into the current tracepoint.
                if len > max_diff {
                    // If the indel is too long, flush any pending segment first.
                    if cur_a_len > 0 || cur_b_len > 0 {
                        tracepoints.push((cur_a_len, cur_b_len));
                        cur_a_len = 0;
                        cur_b_len = 0;
                        cur_diff = 0;
                    }
                    if op == 'I' {
                        tracepoints.push((len, 0));
                    } else {
                        // op == 'D'
                        tracepoints.push((0, len));
                    }
                } else {
                    // If adding this indel would push the diff over the threshold, flush first.
                    if cur_diff + len > max_diff {
                        tracepoints.push((cur_a_len, cur_b_len));
                        cur_a_len = 0;
                        cur_b_len = 0;
                        cur_diff = 0;
                    }
                    // Then accumulate the entire indel.
                    if op == 'I' {
                        cur_a_len += len;
                    } else {
                        // op == 'D'
                        cur_b_len += len;
                    }
                    cur_diff += len;
                }
            },
            '=' | 'M' => {
                // For match-type ops, simply accumulate since they don't add to diff.
                cur_a_len += len;
                cur_b_len += len;
            },
            _ => {
                // Fallback: error
                panic!("Invalid CIGAR operation: {}", op);
            }
        }
    }
    // Flush any remaining segment.
    if cur_a_len > 0 || cur_b_len > 0 {
        tracepoints.push((cur_a_len, cur_b_len));
    }
    tracepoints
}

/// Convert a CIGAR string into tracepoints with diagonal tracking.
/// 
/// Similar to cigar_to_tracepoints but adds diagonal boundary tracking.
/// The diagonal position tracks the relative offset between sequences and helps
/// optimize subsequent alignment by constraining the search space.
/// 
/// @param cigar: The CIGAR string to process
/// @param max_diff: Maximum number of differences allowed in each segment
/// @return Vector of tracepoints with diagonal boundaries: (a_len, b_len, (min_diagonal, max_diagonal))
pub fn cigar_to_banded_tracepoints(
    cigar: &str,
    max_diff: usize,
) -> Vec<(usize, usize, (isize, isize))> {
    let ops = cigar_str_to_cigar_ops(cigar);
    let mut tracepoints = Vec::new();
    // current segment counters:
    let mut cur_a_len = 0;
    let mut cur_b_len = 0;
    let mut cur_diff = 0;

    let mut current_diagonal : isize = 0;  // Current diagonal position
    let mut min_diagonal : isize = 0;      // Lowest diagonal reached
    let mut max_diagonal : isize = 0;      // Highest diagonal reached

    for (mut len, op) in ops {
        match op {
            'X' => {
                // X is splittable
                while len > 0 {
                    let remaining = max_diff.saturating_sub(cur_diff);
                    let step = min(len, remaining);
                    cur_a_len += step;
                    cur_b_len += step;
                    cur_diff += step;
                    len -= step;
                    if cur_diff == max_diff {
                        tracepoints.push((cur_a_len, cur_b_len, (min_diagonal, max_diagonal)));
                        cur_a_len = 0;
                        cur_b_len = 0;
                        cur_diff = 0;
                        current_diagonal = 0;
                        min_diagonal = 0;
                        max_diagonal = 0;
                    }
                }
            },
            'I' | 'D' => {
                // For indels, which are unsplittable, try to incorporate into the current tracepoint.
                if len > max_diff {
                    // If the indel is too long, flush any pending segment first.
                    if cur_a_len > 0 || cur_b_len > 0 {
                        tracepoints.push((cur_a_len, cur_b_len, (min_diagonal, max_diagonal)));
                        cur_a_len = 0;
                        cur_b_len = 0;
                        cur_diff = 0;
                        current_diagonal = 0;
                        min_diagonal = 0;
                        max_diagonal = 0;
                    }
                    if op == 'I' {
                        tracepoints.push((len, 0, (0, len as isize)));
                    } else {
                        // op == 'D'
                        tracepoints.push((0, len, (len as isize, 0)));
                    }
                } else {
                    // If adding this indel would push the diff over the threshold, flush first.
                    if cur_diff + len > max_diff {
                        tracepoints.push((cur_a_len, cur_b_len, (min_diagonal, max_diagonal)));
                        cur_a_len = 0;
                        cur_b_len = 0;
                        cur_diff = 0;
                        current_diagonal = 0;
                        min_diagonal = 0;
                        max_diagonal = 0;
                    }
                    // Then accumulate the entire indel.
                    if op == 'I' {
                        cur_a_len += len;
                    } else {
                        // op == 'D'
                        cur_b_len += len;
                    }
                    cur_diff += len;

                    if op == 'I' {
                        current_diagonal += len as isize;
                        max_diagonal = max_diagonal.max(current_diagonal);
                    } else {
                        // op == 'D'
                        current_diagonal -= len as isize;
                        min_diagonal = min_diagonal.min(current_diagonal);
                    }
                }
            },
            '=' | 'M' => {
                // For match-type ops, simply accumulate since they don't add to diff.
                cur_a_len += len;
                cur_b_len += len;
            },
            _ => {
                // Fallback: error
                panic!("Invalid CIGAR operation: {}", op);
            }
        }
    }
    // Flush any remaining segment.
    if cur_a_len > 0 || cur_b_len > 0 {
        tracepoints.push((cur_a_len, cur_b_len, (min_diagonal, max_diagonal)));
    }
    tracepoints
}

/// Reconstruct a CIGAR string from tracepoints.
/// 
/// For each tracepoint pair, this function performs a detailed alignment of the corresponding
/// sequence segments using WFA alignment. Special cases are handled efficiently:
/// - Pure insertions (a_len > 0, b_len = 0) are directly converted to 'I' operations
/// - Pure deletions (a_len = 0, b_len > 0) are directly converted to 'D' operations
/// - Mixed segments are realigned using the WFA algorithm
/// 
/// @param tracepoints: Vector of (a_len, b_len) pairs defining segment boundaries
/// @param a_seq: Reference sequence string
/// @param b_seq: Query sequence string
/// @param a_start: Starting position in reference sequence
/// @param b_start: Starting position in query sequence
/// @param mismatch: Penalty for mismatches
/// @param gap_open1: Penalty for opening a gap (first gap type)
/// @param gap_ext1: Penalty for extending a gap (first gap type)
/// @param gap_open2: Penalty for opening a gap (second gap type)
/// @param gap_ext2: Penalty for extending a gap (second gap type)
/// @return Reconstructed CIGAR string
pub fn tracepoints_to_cigar(
    tracepoints: &[(usize, usize)],
    a_seq: &[u8],
    b_seq: &[u8],
    a_start: usize,
    b_start: usize,
    mismatch: i32,
    gap_open1: i32,
    gap_ext1: i32,
    gap_open2: i32,
    gap_ext2: i32,
) -> String {
    // Create aligner and configure settings
    let mut aligner = AffineWavefronts::with_penalties_affine2p(0, mismatch, gap_open1, gap_ext1, gap_open2, gap_ext2);

    let mut cigar_ops = Vec::new();
    let mut current_a = a_start;
    let mut current_b = b_start;

    for &(a_len, b_len) in tracepoints {
        // Special case: long indel.
        if a_len > 0 && b_len == 0 {
            // This is an insertion.
            cigar_ops.push((a_len, 'I'));
            current_a += a_len;
        } else if b_len > 0 && a_len == 0 {
            // This is a deletion.
            cigar_ops.push((b_len, 'D'));
            current_b += b_len;
        } else {
            let a_end = current_a + a_len;
            let b_end = current_b + b_len;
            let seg_ops = align_sequences_wfa(
                &a_seq[current_a..a_end],
                &b_seq[current_b..b_end],
                &mut aligner
            );
            cigar_ops.extend(seg_ops);
            current_a = a_end;
            current_b = b_end;
        }
    }
    let merged = merge_cigar_ops(cigar_ops);
    cigar_ops_to_cigar_string(&merged)
}

/// Reconstruct a CIGAR string from banded tracepoints.
/// 
/// Similar to tracepoints_to_cigar but utilizes the diagonal boundary
/// information to constrain the WFA alignment search space. This improves performance
/// by limiting the search to the banded region where the alignment is expected to be found.
/// 
/// @param tracepoints: Vector of (a_len, b_len, (min_diagonal, max_diagonal)) triples
/// @param a_seq: Reference sequence string
/// @param b_seq: Query sequence string
/// @param a_start: Starting position in reference sequence
/// @param b_start: Starting position in query sequence
/// @param mismatch: Penalty for mismatches
/// @param gap_open1: Penalty for opening a gap (first gap type)
/// @param gap_ext1: Penalty for extending a gap (first gap type)
/// @param gap_open2: Penalty for opening a gap (second gap type)
/// @param gap_ext2: Penalty for extending a gap (second gap type)
/// @return Reconstructed CIGAR string
pub fn banded_tracepoints_to_cigar(
    tracepoints: &[(usize, usize, (isize, isize))],
    a_seq: &[u8],
    b_seq: &[u8],
    a_start: usize,
    b_start: usize,
    mismatch: i32,
    gap_open1: i32,
    gap_ext1: i32,
    gap_open2: i32,
    gap_ext2: i32,
) -> String {
    // Create aligner and configure settings
    let mut aligner = AffineWavefronts::with_penalties_affine2p(0, mismatch, gap_open1, gap_ext1, gap_open2, gap_ext2);

    let mut cigar_ops = Vec::new();
    let mut current_a = a_start;
    let mut current_b = b_start;

    for &(a_len, b_len, (min_k, max_k)) in tracepoints {
        // Special case: long indel.
        if a_len > 0 && b_len == 0 {
            // This is an insertion.
            cigar_ops.push((a_len, 'I'));
            current_a += a_len;
        } else if b_len > 0 && a_len == 0 {
            // This is a deletion.
            cigar_ops.push((b_len, 'D'));
            current_b += b_len;
        } else {
            let a_end = current_a + a_len;
            let b_end = current_b + b_len;
            aligner.set_heuristic(&HeuristicStrategy::BandedStatic { band_min_k: (min_k - 1) as i32, band_max_k: (max_k + 1) as i32 });
            let seg_ops = align_sequences_wfa(
                &a_seq[current_a..a_end],
                &b_seq[current_b..b_end],
                &mut aligner
            );
            cigar_ops.extend(seg_ops);
            current_a = a_end;
            current_b = b_end;
        }
    }
    let merged = merge_cigar_ops(cigar_ops);
    cigar_ops_to_cigar_string(&merged)
}


// Helper functions

/// Merge consecutive CIGAR operations of the same type.
/// 
/// @param ops: Vector of (length, operation) pairs
/// @return Vector with adjacent identical operations combined
fn merge_cigar_ops(ops: Vec<(usize, char)>) -> Vec<(usize, char)> {
    if ops.is_empty() {
        return ops;
    }
    let mut merged = Vec::new();
    let (mut count, mut op) = ops[0];
    for &(c, o) in ops.iter().skip(1) {
        if o == op {
            count += c;
        } else {
            merged.push((count, op));
            op = o;
            count = c;
        }
    }
    merged.push((count, op));
    merged
}

/// Parse a CIGAR string into a vector of (length, operation) pairs.
/// 
/// @param cigar: CIGAR string (e.g., "10M2D5M")
/// @return Vector of (length, operation) pairs
fn cigar_str_to_cigar_ops(cigar: &str) -> Vec<(usize, char)> {
    let mut ops = Vec::new();
    let mut num = String::new();
    for ch in cigar.chars() {
        if ch.is_digit(10) {
            num.push(ch);
        } else {
            if let Ok(n) = num.parse::<usize>() {
                ops.push((n, ch));
            }
            num.clear();
        }
    }
    ops
}

/// Convert a byte array of CIGAR operations into a vector of (length, operation) pairs.
/// Handles the special case of 'M' (77 in ASCII) being treated as '='.
/// 
/// @param ops: Byte array of CIGAR operations
/// @return Vector of (length, operation) pairs
fn cigar_u8_to_cigar_ops(ops: &[u8]) -> Vec<(usize, char)> {
    let mut result = Vec::new();
    let mut count = 1;
    let mut current_op = if ops[0] == 77 { '=' } else { ops[0] as char };
    
    for &byte in ops.iter().skip(1) {
        let op = if byte == 77 { '=' } else { byte as char };
        if op == current_op {
            count += 1;
        } else {
            result.push((count, current_op));
            current_op = op;
            count = 1;
        }
    }
    
    // Push the final operation
    result.push((count, current_op));
    
    result
}

/// Format a vector of CIGAR operations as a standard CIGAR string.
/// 
/// @param ops: Vector of (length, operation) pairs
/// @return Formatted CIGAR string (e.g., "10M2D5M")
fn cigar_ops_to_cigar_string(ops: &[(usize, char)]) -> String {
    ops.iter()
        .map(|(len, op)| format!("{}{}", len, op))
        .collect::<Vec<_>>()
        .join("")
}

/// Align two sequence segments using WFA algorithm.
/// 
/// @param a: Reference sequence segment
/// @param b: Query sequence segment
/// @param aligner: Configured WFA aligner
/// @return Vector of CIGAR operations for the alignment
fn align_sequences_wfa(
    a: &[u8],
    b: &[u8],
    aligner: &mut AffineWavefronts
) -> Vec<(usize, char)> {   
    // Do the alignment b vs a (it is not a typo) to have insertions/deletions in the query as Is/Ds in the CIGAR string
    let status = aligner.align(b, a);
    
    match status {
        AlignmentStatus::Completed => {
            cigar_u8_to_cigar_ops(aligner.cigar())
        },
        s => {
            panic!("Alignment failed with status: {:?}", s);
        }
    }
}

// /// With the inverted logic, a is consumed by insertions:
// /// so =, X, and I (and M) consume A.
// fn consumes_a(op: char) -> bool {
//     op == '=' || op == 'X' || op == 'I' || op == 'M'
// }

// /// With the inverted logic, b is consumed by deletions:
// /// so =, X, and D (and M) consume B.
// fn consumes_b(op: char) -> bool {
//     op == '=' || op == 'X' || op == 'D' || op == 'M'
// }

// /// Convert a CIGAR string into tracepoints.
// /// Given:
// /// - `cigar`: the alignment CIGAR string in extended format.
// /// - `a_start`, `a_end`: the positions on sequence A that the alignment spans.
// /// - `b_start`, `b_end`: similarly for sequence B.
// /// - `delta`: the spacing for tracepoints.
// /// Returns a vector of (diff_count, b_bases) pairs for each A–segment.
// /// BUGGY!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// fn cigar_to_tracepoints(
//     cigar: &str,
//     a_start: usize,
//     a_end: usize,
//     delta: usize,
// ) -> Vec<(usize, usize)> {
//     let ops = cigar_str_to_cigar_ops(cigar);

//     // Calculate the next A boundary (tracepoint) after a_start.
//     let mut next_thresh = if a_start % delta == 0 {
//         a_start + delta
//     } else {
//         ((a_start / delta) + 1) * delta
//     };

//     let mut a_pos = a_start;
//     let mut current_diffs = 0;
//     let mut current_b_bases = 0;

//     let mut tracepoints = Vec::new();

//     for (mut len, op) in ops {
//         // For op that does not consume A (e.g. deletions) we process it wholly.
//         if !consumes_a(op) {
//             if consumes_b(op) {
//                 current_b_bases += len;
//             }
//             if is_edit(op) {
//                 current_diffs += len;
//             }
//             continue;
//         }
//         // For op that consumes A, we may need to “split” it if it spans a tracepoint boundary.
//         while len > 0 {
//             let remaining_to_thresh = next_thresh.saturating_sub(a_pos);
//             let step = min(len, remaining_to_thresh);

//             // Determine how many bases on B this step consumes.
//             let b_step = if consumes_b(op) { step } else { 0 };

//             // Update counters.
//             a_pos += step;
//             current_b_bases += b_step;
//             if is_edit(op) {
//                 current_diffs += step;
//             }

//             len -= step;

//             // If we’ve reached a tracepoint boundary (i.e. a_pos == next_thresh)
//             if a_pos == next_thresh {
//                 // Record the tracepoint for this segment.
//                 tracepoints.push((current_diffs, current_b_bases));
//                 // Reset accumulators for the next segment.
//                 current_diffs = 0;
//                 current_b_bases = 0;
//                 next_thresh += delta;
//             }
//         }
//     }
//     // Always record a final segment if we haven't reached a_end or if there is leftover.
//     if a_pos < a_end || current_diffs != 0 || current_b_bases != 0 {
//         tracepoints.push((current_diffs, current_b_bases));
//     }
//     tracepoints
// }

// /// Given tracepoints (the per-segment (d,b) pairs), re‐construct the full CIGAR string.
// /// To “fill in” the alignment for each segment we use a global alignment with affine gap penalties
// /// on the corresponding sub–sequences. The A boundaries are determined by a_start, a_end and delta;
// /// for each segment the B–interval length is taken from the tracepoint b value.
// ///
// /// Note: here we ignore the d value (edit count) from the tracepoint.
// fn tracepoints_to_cigar(
//     tracepoints: &[(usize, usize)],
//     a_seq: &str,
//     b_seq: &str,
//     a_start: usize,
//     a_end: usize,
//     b_start: usize,
//     _b_end: usize,
//     delta: usize,
//     mismatch: i32,
//     gap_open_i: i32,
//     gap_extend_i: i32,
//     gap_open_d: i32,
//     gap_extend_d: i32,
// ) -> String {
//     // First, compute the number of intervals based solely on A consumption.
//     // (This does not count any trailing B-only operations.)
//     let consumed_intervals = ((a_end - a_start) + delta - 1) / delta;
//     // If we recorded more tracepoints than that, we have a trailing segment.
//     let extra = if tracepoints.len() > consumed_intervals { 1 } else { 0 };
//     //let total_intervals = consumed_intervals + extra;

//     // Now compute the boundaries.
//     // For the consumed intervals, we use the standard boundaries.
//     let mut boundaries = Vec::new();
//     boundaries.push(a_start);
//     let mut next = if a_start % delta == 0 {
//         a_start + delta
//     } else {
//         ((a_start / delta) + 1) * delta
//     };
//     while boundaries.len() < consumed_intervals {
//         boundaries.push(next);
//         next += delta;
//     }
//     boundaries.push(a_end);
//     // If there's an extra tracepoint (trailing B-only segment),
//     // add an extra boundary equal to a_end so that the final interval is zero-length in A.
//     if extra == 1 {
//         boundaries.push(a_end);
//     }
//     // Now, boundaries.len() - 1 should equal total_intervals (which equals tracepoints.len()).
//     if boundaries.len() - 1 != tracepoints.len() {
//         panic!(
//             "Mismatch: {} intervals vs {} tracepoint pairs",
//             boundaries.len() - 1,
//             tracepoints.len()
//         );
//     }

//     // Create aligner and configure settings
//     let mut aligner = AffineWavefronts::with_penalties_affine2p(0, mismatch, gap_open_i, gap_extend_i, gap_open_d, gap_extend_d);

//     let mut cigar_ops = Vec::new();
//     let mut current_b = b_start;
//     for (i, &(_d, b_len)) in tracepoints.iter().enumerate() {
//         let a_left = boundaries[i];
//         let a_right = boundaries[i + 1];

//         // Extract the sub–sequences.
//         let a_sub = &a_seq[a_left..a_right];
//         let b_sub = &b_seq[current_b..current_b + (b_len as usize)];

//         // Add debug output
//         // eprintln!("Segment {}: A[{}..{}] (len {}), B[{}..{}] (len {})",
//         // i, a_left, a_right, a_right - a_left,
//         // current_b, current_b + b_len, b_len);

//         // Align the two segments using affine gap penalties.
//         let mut seg_cigar = align_sequences_wfa(
//             a_sub,
//             b_sub,
//             &mut aligner
//         );

//         // Append to our overall CIGAR operations.
//         cigar_ops.append(&mut seg_cigar);
//         current_b += b_len as usize;
//     }

//     // Merge adjacent operations if needed.
//     let merged = merge_cigar_ops(cigar_ops);
//     cigar_ops_to_cigar_string(&merged)
// }

// /// Returns true if the op counts as an edit (difference).
// fn is_edit(op: char) -> bool {
//     op == 'X' || op == 'I' || op == 'D'
// }

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tracepoint_generation() {
        // Test CIGAR string
        let cigar = "10M2D5M2I3M";
        
        // Define test cases: (max_diff, expected_tracepoints, expected_banded_tracepoints)
        let test_cases = vec![
            // Case 1: No differences allowed - each operation becomes its own segment
            (0, 
             vec![(10, 10), (0, 2), (5, 5), (2, 0), (3, 3)],
             vec![(10, 10, (0, 0)), (0, 2, (2, 0)), (5, 5, (0, 0)), (2, 0, (0, 2)), (3, 3, (0, 0))]),
            
            // Case 2: Allow up to 2 differences in each segment
            (2,
             vec![(15, 17), (5, 3)],
             vec![(15, 17, (-2, 0)), (5, 3, (0, 2))]),
            
            // Case 3: Allow up to 5 differences - combines all operations
            (5,
             vec![(20, 20)],
             vec![(20, 20, (-2, 0))])
        ];
        
        // Run each test case
        for (i, (max_diff, expected_tracepoints, expected_banded_tracepoints)) in test_cases.iter().enumerate() {
            // Get actual results
            let tracepoints = cigar_to_tracepoints(&cigar, *max_diff);
            let banded_tracepoints = cigar_to_banded_tracepoints(&cigar, *max_diff);
            
            // Check basic tracepoints
            assert_eq!(tracepoints, *expected_tracepoints,
                       "Test case {}: Basic tracepoints with max_diff={} incorrect", i+1, max_diff);
                       
            // Check banded tracepoints
            assert_eq!(banded_tracepoints, *expected_banded_tracepoints,
                       "Test case {}: Banded tracepoints with max_diff={} incorrect", i+1, max_diff);
            
            // Verify both implementations are consistent in terms of segment lengths
            assert_eq!(tracepoints.len(), banded_tracepoints.len(), 
                       "Test case {}: Implementations should produce the same number of segments", i+1);
            
            for (j, ((a_len, b_len), (a_len_banded, b_len_banded, _))) in 
                tracepoints.iter().zip(banded_tracepoints.iter()).enumerate() {
                assert_eq!(
                    (a_len, b_len), 
                    (a_len_banded, b_len_banded),
                    "Test case {}, segment {}: Length mismatch - Standard: ({}, {}) vs Banded: ({}, {})",
                    i+1, j, a_len, b_len, a_len_banded, b_len_banded
                );
            }
        }
    }

    #[test]
    fn test_cigar_roundtrip() {
        let original_cigar = "1=1I18=";
        let a_seq = b"ACGTACGTACACGTACGTAC";  // 20 bases
        let b_seq = b"AGTACGTACACGTACGTAC";   // 19 bases (missing C)
        let max_diff = 5;
        
        // Test basic tracepoints
        let tracepoints = cigar_to_tracepoints(&original_cigar, max_diff);
        let basic_cigar = tracepoints_to_cigar(
            &tracepoints,
            a_seq, b_seq, 
            0, 0,  // sequences and start positions
            2, 4, 2, 6, 1        // alignment penalties
        );
        assert_eq!(basic_cigar, original_cigar, "Basic implementation failed");
        
        // Test banded tracepoints
        let banded_tracepoints: Vec<(usize, usize, (isize, isize))> = cigar_to_banded_tracepoints(&original_cigar, max_diff);
        let banded_cigar = banded_tracepoints_to_cigar(
            &banded_tracepoints,
            a_seq, b_seq, 
            0, 0,  // sequences and start positions
            2, 4, 2, 6, 1        // alignment penalties
        );
        assert_eq!(banded_cigar, original_cigar, "Banded implementation failed");
    }
}
