use std::cmp::min;
use lib_wfa2::affine_wavefront::{AffineWavefronts, AlignmentStatus, HeuristicStrategy};

/// Convert a CIGAR string into variable–delta tracepoints.
/// Instead of using fixed A–intervals, we accumulate bases (and differences)
/// until we reach a diff threshold of max_diff. For match-like ops ('=', 'M', and 'X'),
/// we split as needed so that the diff count never exceeds max_diff.
/// For indels ('I' and 'D'), we incorporate them into the current tracepoint
/// if they are short enough. If adding the indel would exceed the threshold,
/// we flush the current segment. Indels are unsplittable.
pub fn cigar_to_tracepoints_variable(
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

pub fn cigar_to_tracepoints_variable2(
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


/// Reconstruct a CIGAR string from variable-delta tracepoints.
/// For each tracepoint, we simply generate an insertion or deletion op by
/// inspecting the saved A and B bases. /// Otherwise, we realign the
/// corresponding segments using either heuristic or full WFA alignment.
pub fn tracepoints_to_cigar_variable(
    tracepoints: &[(usize, usize)],
    a_seq: &str,
    b_seq: &str,
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

pub fn tracepoints_to_cigar_variable2(
    tracepoints: &[(usize, usize, (isize, isize))],
    a_seq: &str,
    b_seq: &str,
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

/// With the inverted logic, a is consumed by insertions:
/// so =, X, and I (and M) consume A.
fn consumes_a(op: char) -> bool {
    op == '=' || op == 'X' || op == 'I' || op == 'M'
}

/// With the inverted logic, b is consumed by deletions:
/// so =, X, and D (and M) consume B.
fn consumes_b(op: char) -> bool {
    op == '=' || op == 'X' || op == 'D' || op == 'M'
}

/// Helper: Merge two vectors of CIGAR operations (merging adjacent ops of the same kind).
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

/// Parse a CIGAR string into a vector of (length, op) pairs.
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

/// Convert a vector of (length, op) pairs to a CIGAR string.
fn cigar_ops_to_cigar_string(ops: &[(usize, char)]) -> String {
    ops.iter()
        .map(|(len, op)| format!("{}{}", len, op))
        .collect::<Vec<_>>()
        .join("")
}

fn align_sequences_wfa(
    a: &str, 
    b: &str,
    aligner: &mut AffineWavefronts
) -> Vec<(usize, char)> {   
    // Do the alignment b vs a (it is not a typo) to have insertions/deletions in the query as Is/Ds in the CIGAR string
    let status = aligner.align(b.as_bytes(), a.as_bytes());
    
    match status {
        AlignmentStatus::Completed => {
            cigar_u8_to_cigar_ops(aligner.cigar())
        },
        s => {
            panic!("Alignment failed with status: {:?}", s);
        }
    }
}

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
    fn it_works() {
        let cigar = "10M2D5M2I3M";

        let tracepoints = cigar_to_tracepoints_variable(&cigar, 0);
        assert_eq!(tracepoints, vec![(10, 10), (0, 2), (5, 5), (2, 0), (3, 3)]);

        let tracepoints = cigar_to_tracepoints_variable(&cigar, 2);
        assert_eq!(tracepoints, vec![(15, 17), (5, 3)]);

        let tracepoints = cigar_to_tracepoints_variable(&cigar, 5);
        assert_eq!(tracepoints, vec![(20, 20)]);
    }

    #[test]
    fn xxx() {
        let cigar = "1=1I18=";
        let tracepoints = cigar_to_tracepoints_variable(&cigar, 5);
        
        // Reconstruct CIGAR from tracepoints
        let a_seq = "ACGTACGTACACGTACGTAC";
        let b_seq = "AGTACGTACACGTACGTAC";
        let reconstructed_cigar = tracepoints_to_cigar_variable(
            &tracepoints,
            a_seq,
            b_seq,
            0,  // a_start
            0,  // b_start
            2,  // mismatch penalty
            4,  // gap_open1
            2,  // gap_ext1
            6,  // gap_open2
            1,  // gap_ext2
        );

        assert_eq!(reconstructed_cigar, cigar);
    }
}
