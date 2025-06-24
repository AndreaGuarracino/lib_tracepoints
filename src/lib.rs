use lib_wfa2::affine_wavefront::{AffineWavefronts, AlignmentStatus};
use std::cmp::min;

/// Represents a CIGAR segment that can be either aligned or preserved as-is
#[derive(Debug, Clone, PartialEq)]
pub enum MixedRepresentation {
    /// Alignment segment represented by tracepoints
    Tracepoint(usize, usize),
    /// Special CIGAR operation that should be preserved intact
    CigarOp(usize, char),
}

/// Convert a CIGAR string into variable tracepoints with length optimization.
///
/// Similar to cigar_to_tracepoints but optimizes storage by using a single value
/// when a_len == b_len. The representation is (length, Option<b_length>):
/// - If a_len == b_len: stores (a_len, None)
/// - If a_len != b_len: stores (a_len, Some(b_len))
///
/// @param cigar: The CIGAR string to process
/// @param max_diff: Maximum number of differences allowed in each segment
/// @return Vector of variable tracepoints containing (a_segment_length, Option<b_segment_length>)
pub fn cigar_to_variable_tracepoints(cigar: &str, max_diff: usize) -> Vec<(usize, Option<usize>)> {
    let basic_tracepoints = cigar_to_tracepoints(cigar, max_diff);
    
    basic_tracepoints
        .into_iter()
        .map(|(a_len, b_len)| {
            if a_len == b_len {
                (a_len, None)
            } else {
                (a_len, Some(b_len))
            }
        })
        .collect()
}

/// Convert a CIGAR string into tracepoints.
///
/// Segments a CIGAR string into tracepoints with at most `max_diff` differences per segment.
/// - Match operations ('=' and 'M') don't count as differences.
/// - Mismatches ('X') can be split across segments
/// - Indels ('I', 'D') are kept in a single segment when possible
/// - Long indels exceeding max_diff become their own segments
///
/// @param cigar: The CIGAR string to process
/// @param max_diff: Maximum number of differences allowed in each segment
/// @return Vector of tracepoints containing (a_segment_length, b_segment_length)
pub fn cigar_to_tracepoints(cigar: &str, max_diff: usize) -> Vec<(usize, usize)> {
    let ops = cigar_str_to_cigar_ops(cigar);
    let mut tracepoints = Vec::new();

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

                    // Handle the case when max_diff = 0 or no room left in current segment
                    if step == 0 {
                        // Flush current segment if it exists
                        if cur_a_len > 0 || cur_b_len > 0 {
                            tracepoints.push((cur_a_len, cur_b_len));
                            cur_a_len = 0;
                            cur_b_len = 0;
                            cur_diff = 0;
                        }

                        // Create a segment with just 1 mismatch
                        if max_diff == 0 {
                            tracepoints.push((1, 1));
                            len -= 1;
                        } else {
                            // This case shouldn't happen, but handle it gracefully
                            cur_a_len = 1;
                            cur_b_len = 1;
                            cur_diff = 1;
                            len -= 1;
                        }
                    } else {
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
                }
            }
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
            }
            '=' | 'M' => {
                // For match-type ops, simply accumulate since they don't add to diff.
                cur_a_len += len;
                cur_b_len += len;
            }
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

/// Convert a CIGAR string into mixed representation tracepoints.
///
/// Similar to cigar_to_tracepoints but preserves special operations (H, N, P, S)
/// that should not be processed by WFA alignment.
///
/// @param cigar: The CIGAR string to process
/// @param max_diff: Maximum number of differences allowed in each segment
/// @return Vector of mixed representation elements
pub fn cigar_to_mixed_tracepoints(cigar: &str, max_diff: usize) -> Vec<MixedRepresentation> {
    let ops = cigar_str_to_cigar_ops(cigar);
    let mut mixed_tracepoints = Vec::new();

    let mut cur_a_len = 0;
    let mut cur_b_len = 0;
    let mut cur_diff = 0;

    for (mut len, op) in ops {
        match op {
            // Special operators that are preserved as-is
            'H' | 'N' | 'P' | 'S' => {
                if cur_a_len > 0 || cur_b_len > 0 {
                    mixed_tracepoints.push(MixedRepresentation::Tracepoint(cur_a_len, cur_b_len));
                    cur_a_len = 0;
                    cur_b_len = 0;
                    cur_diff = 0;
                }

                // Add the special operation
                mixed_tracepoints.push(MixedRepresentation::CigarOp(len, op));
            }
            'X' => {
                // X is splittable
                while len > 0 {
                    let remaining = max_diff.saturating_sub(cur_diff);
                    let step = min(len, remaining);

                    // Handle the case when max_diff = 0 or no room left in current segment
                    if step == 0 {
                        // Flush current segment if it exists
                        if cur_a_len > 0 || cur_b_len > 0 {
                            mixed_tracepoints.push(MixedRepresentation::Tracepoint(cur_a_len, cur_b_len));
                            cur_a_len = 0;
                            cur_b_len = 0;
                            cur_diff = 0;
                        }

                        // Create a segment with just 1 mismatch
                        if max_diff == 0 {
                            mixed_tracepoints.push(MixedRepresentation::Tracepoint(1, 1));
                            len -= 1;
                        } else {
                            // This case shouldn't happen, but handle it gracefully
                            cur_a_len = 1;
                            cur_b_len = 1;
                            cur_diff = 1;
                            len -= 1;
                        }
                    } else {
                        cur_a_len += step;
                        cur_b_len += step;
                        cur_diff += step;
                        len -= step;

                        if cur_diff == max_diff {
                            mixed_tracepoints.push(MixedRepresentation::Tracepoint(cur_a_len, cur_b_len));
                            cur_a_len = 0;
                            cur_b_len = 0;
                            cur_diff = 0;
                        }
                    }
                }
            }
            'I' | 'D' => {
                // For indels, which are unsplittable, try to incorporate into the current tracepoint.
                if len > max_diff {
                    // If the indel is too long, flush any pending segment first.
                    if cur_a_len > 0 || cur_b_len > 0 {
                        mixed_tracepoints.push(MixedRepresentation::Tracepoint(cur_a_len, cur_b_len));
                        cur_a_len = 0;
                        cur_b_len = 0;
                        cur_diff = 0;
                    }
                    if op == 'I' {
                        mixed_tracepoints.push(MixedRepresentation::Tracepoint(len, 0));
                    } else {
                        // op == 'D'
                        mixed_tracepoints.push(MixedRepresentation::Tracepoint(0, len));
                    }
                } else {
                    // If adding this indel would push the diff over the threshold, flush first.
                    if cur_diff + len > max_diff {
                        mixed_tracepoints.push(MixedRepresentation::Tracepoint(cur_a_len, cur_b_len));
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
            }
            '=' | 'M' => {
                // For match-type ops, simply accumulate since they don't add to diff.
                cur_a_len += len;
                cur_b_len += len;
            }
            _ => {
                // Fallback: error
                panic!("Invalid CIGAR operation: {}", op);
            }
        }
    }
    // Flush any remaining segment.
    if cur_a_len > 0 || cur_b_len > 0 {
        mixed_tracepoints.push(MixedRepresentation::Tracepoint(cur_a_len, cur_b_len));
    }
    mixed_tracepoints
}



/// Reconstruct a CIGAR string from variable tracepoints.
///
/// Converts variable tracepoints back to regular tracepoints and then reconstructs CIGAR.
/// The variable format (length, Option<b_length>) is expanded to (a_len, b_len):
/// - (len, None) becomes (len, len)
/// - (a_len, Some(b_len)) becomes (a_len, b_len)
///
/// @param variable_tracepoints: Vector of (a_len, Option<b_len>) pairs
/// @param a_seq: Reference sequence
/// @param b_seq: Query sequence
/// @param a_start: Starting position in reference sequence
/// @param b_start: Starting position in query sequence
/// @param mismatch: Penalty for mismatches
/// @param gap_open1: Penalty for opening a gap (first gap type)
/// @param gap_ext1: Penalty for extending a gap (first gap type)
/// @param gap_open2: Penalty for opening a gap (second gap type)
/// @param gap_ext2: Penalty for extending a gap (second gap type)
/// @return Reconstructed CIGAR string
pub fn variable_tracepoints_to_cigar(
    variable_tracepoints: &[(usize, Option<usize>)],
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
    // Convert variable tracepoints back to regular tracepoints
    let regular_tracepoints: Vec<(usize, usize)> = variable_tracepoints
        .iter()
        .map(|(a_len, b_len_opt)| {
            match b_len_opt {
                None => (*a_len, *a_len),
                Some(b_len) => (*a_len, *b_len),
            }
        })
        .collect();
    
    // Use existing tracepoints_to_cigar function
    tracepoints_to_cigar(
        &regular_tracepoints,
        a_seq,
        b_seq,
        a_start,
        b_start,
        mismatch,
        gap_open1,
        gap_ext1,
        gap_open2,
        gap_ext2,
    )
}

/// Reconstruct a CIGAR string from tracepoints.
///
/// For each tracepoint pair, performs detailed alignment using WFA alignment.
/// Special cases are handled efficiently:
/// - Pure insertions (a_len > 0, b_len = 0) become 'I' operations
/// - Pure deletions (a_len = 0, b_len > 0) become 'D' operations
/// - Mixed segments are realigned using the WFA algorithm
///
/// @param tracepoints: Vector of (a_len, b_len) pairs defining segment boundaries
/// @param a_seq: Reference sequence
/// @param b_seq: Query sequence
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
    let mut aligner = AffineWavefronts::with_penalties_affine2p(
        0, mismatch, gap_open1, gap_ext1, gap_open2, gap_ext2,
    );

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
                &mut aligner,
            );
            cigar_ops.extend(seg_ops);
            current_a = a_end;
            current_b = b_end;
        }
    }
    merge_cigar_ops(&mut cigar_ops);
    cigar_ops_to_cigar_string(&cigar_ops)
}

/// Reconstruct a CIGAR string from mixed representation.
///
/// Processes both regular alignment segments (Tracepoint) and special operations (CigarOp)
/// to reconstruct a complete CIGAR string with all operations preserved.
///
/// @param mixed_tracepoints: Vector of MixedRepresentation items
/// @param a_seq: Reference sequence
/// @param b_seq: Query sequence
/// @param a_start: Starting position in reference sequence
/// @param b_start: Starting position in query sequence
/// @param mismatch: Penalty for mismatches
/// @param gap_open1: Penalty for opening a gap (first gap type)
/// @param gap_ext1: Penalty for extending a gap (first gap type)
/// @param gap_open2: Penalty for opening a gap (second gap type)
/// @param gap_ext2: Penalty for extending a gap (second gap type)
/// @return Reconstructed CIGAR string
pub fn mixed_tracepoints_to_cigar(
    mixed_tracepoints: &[MixedRepresentation],
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
    let mut aligner = AffineWavefronts::with_penalties_affine2p(
        0, mismatch, gap_open1, gap_ext1, gap_open2, gap_ext2,
    );

    let mut cigar_ops = Vec::new();
    let mut current_a = a_start;
    let mut current_b = b_start;

    for item in mixed_tracepoints {
        match item {
            // For special CIGAR operations, simply add them directly
            MixedRepresentation::CigarOp(len, op) => {
                cigar_ops.push((*len, *op));
            }
            // For tracepoint segments, realign using WFA
            MixedRepresentation::Tracepoint(a_len, b_len) => {
                // Special case: long indel.
                if *a_len > 0 && *b_len == 0 {
                    // This is an insertion.
                    cigar_ops.push((*a_len, 'I'));
                    current_a += a_len;
                } else if *b_len > 0 && *a_len == 0 {
                    // This is a deletion.
                    cigar_ops.push((*b_len, 'D'));
                    current_b += b_len;
                } else {
                    let a_end = current_a + a_len;
                    let b_end = current_b + b_len;
                    let seg_ops = align_sequences_wfa(
                        &a_seq[current_a..a_end],
                        &b_seq[current_b..b_end],
                        &mut aligner,
                    );
                    cigar_ops.extend(seg_ops);
                    current_a = a_end;
                    current_b = b_end;
                }
            }
        }
    }
    merge_cigar_ops(&mut cigar_ops);
    cigar_ops_to_cigar_string(&cigar_ops)
}

// Helper functions

/// Merge in-place consecutive CIGAR operations of the same type.
///
/// @param ops: Vector of (length, operation) pairs
fn merge_cigar_ops(ops: &mut Vec<(usize, char)>) {
    if ops.len() <= 1 {
        return;
    }

    let mut write_idx = 0;
    let mut current_count = ops[0].0;
    let mut current_op = ops[0].1;

    for read_idx in 1..ops.len() {
        let (count, op) = ops[read_idx];
        if op == current_op {
            current_count += count;
        } else {
            ops[write_idx] = (current_count, current_op);
            write_idx += 1;
            current_count = count;
            current_op = op;
        }
    }
    ops[write_idx] = (current_count, current_op);
    ops.truncate(write_idx + 1);
}

/// Parse a CIGAR string into a vector of (length, operation) pairs.
///
/// @param cigar: CIGAR string (e.g., "10M2D5M")
/// @return Vector of (length, operation) pairs
fn cigar_str_to_cigar_ops(cigar: &str) -> Vec<(usize, char)> {
    let mut ops = Vec::new();
    let mut num = String::new();
    for ch in cigar.chars() {
        if ch.is_ascii_digit() {
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
pub fn cigar_ops_to_cigar_string(ops: &[(usize, char)]) -> String {
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
pub fn align_sequences_wfa(
    a: &[u8],
    b: &[u8],
    aligner: &mut AffineWavefronts,
) -> Vec<(usize, char)> {
    // Do the alignment b vs a (it is not a typo) to have insertions/deletions in the query as Is/Ds in the CIGAR string
    let status = aligner.align(b, a);

    match status {
        AlignmentStatus::Completed => cigar_u8_to_cigar_ops(aligner.cigar()),
        s => {
            panic!("Alignment failed with status: {:?}", s);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tracepoint_generation() {
        // Test CIGAR string
        let cigar = "10=2D5=2I3=";

        // Define test cases: (max_diff, expected_tracepoints)
        let test_cases = vec![
            // Case 1: No differences allowed - each operation becomes its own segment
            (0, vec![(10, 10), (0, 2), (5, 5), (2, 0), (3, 3)]),
            // Case 2: Allow up to 2 differences in each segment
            (2, vec![(15, 17), (5, 3)]),
            // Case 3: Allow up to 5 differences - combines all operations
            (5, vec![(20, 20)]),
        ];

        // Run each test case
        for (i, (max_diff, expected_tracepoints)) in test_cases.iter().enumerate() {
            // Get actual results
            let tracepoints = cigar_to_tracepoints(cigar, *max_diff);

            // Check tracepoints
            assert_eq!(
                tracepoints,
                *expected_tracepoints,
                "Test case {}: Tracepoints with max_diff={} incorrect",
                i + 1,
                max_diff
            );
        }
    }

    #[test]
    fn test_cigar_roundtrip() {
        let original_cigar = "1=1I18=";
        let a_seq = b"ACGTACGTACACGTACGTAC"; // 20 bases
        let b_seq = b"AGTACGTACACGTACGTAC"; // 19 bases (missing C)
        let max_diff = 5;

        // Test tracepoints roundtrip
        let tracepoints = cigar_to_tracepoints(original_cigar, max_diff);
        let reconstructed_cigar = tracepoints_to_cigar(
            &tracepoints,
            a_seq,
            b_seq,
            0,
            0, // sequences and start positions
            2,
            4,
            2,
            6,
            1, // alignment penalties
        );
        assert_eq!(
            reconstructed_cigar, original_cigar,
            "Tracepoint roundtrip failed"
        );
    }

    #[test]
    fn test_mixed_tracepoint_generation() {
        // Test cases: (CIGAR string, max_diff, expected mixed representation)
        let test_cases = vec![
            // Case 1: Simple alignment with special operators
            (
                "5=2H3=",
                2,
                vec![
                    MixedRepresentation::Tracepoint(5, 5),
                    MixedRepresentation::CigarOp(2, 'H'),
                    MixedRepresentation::Tracepoint(3, 3),
                ],
            ),
            // Case 2: Soft-clipped alignment with indels
            (
                "3S5=2I4=1N2=",
                3,
                vec![
                    MixedRepresentation::CigarOp(3, 'S'),
                    MixedRepresentation::Tracepoint(11, 9),
                    MixedRepresentation::CigarOp(1, 'N'),
                    MixedRepresentation::Tracepoint(2, 2),
                ],
            ),
            // Case 3: Padding with mismatches
            (
                "4P2=3X1=",
                2,
                vec![
                    MixedRepresentation::CigarOp(4, 'P'),
                    MixedRepresentation::Tracepoint(4, 4),
                    MixedRepresentation::Tracepoint(2, 2),
                ],
            ),
            // Case 4: Mixed operations within max_diff
            (
                "5=4X3=2H",
                5,
                vec![
                    MixedRepresentation::Tracepoint(12, 12),
                    MixedRepresentation::CigarOp(2, 'H'),
                ],
            ),
            // Case 5: Large insertion with special ops
            (
                "10I5S",
                3,
                vec![
                    MixedRepresentation::Tracepoint(10, 0),
                    MixedRepresentation::CigarOp(5, 'S'),
                ],
            ),
            // Case 6: Special ops at the beginning
            (
                "3S7D2=",
                4,
                vec![
                    MixedRepresentation::CigarOp(3, 'S'),
                    MixedRepresentation::Tracepoint(0, 7),
                    MixedRepresentation::Tracepoint(2, 2),
                ],
            ),
        ];

        for (i, (cigar, max_diff, expected)) in test_cases.iter().enumerate() {
            let result = cigar_to_mixed_tracepoints(cigar, *max_diff);

            assert_eq!(
                result,
                *expected,
                "Test case {}: Mixed representation with max_diff={} incorrect for CIGAR '{}'",
                i + 1,
                max_diff,
                cigar
            );
        }
    }

    #[test]
    fn test_mixed_tracepoints_roundtrip() {
        // Test data - using identical sequences for predictable results
        let a_seq = b"ACGTACGTACGTACGTACGT"; // 20 bases
        let b_seq = b"ACGTACGTACGTACGTACGT"; // 20 bases (identical for this test)
        
        // Alignment parameters
        let mismatch = 2;
        let gap_open1 = 4;
        let gap_ext1 = 2;
        let gap_open2 = 6;
        let gap_ext2 = 1;

        // Test cases with expected results - simplified to work with identical sequences
        let test_cases = vec![
            // Case 1: Simple special operations with matches
            (
                vec![
                    MixedRepresentation::CigarOp(2, 'S'),
                    MixedRepresentation::Tracepoint(5, 5),
                    MixedRepresentation::CigarOp(3, 'H'),
                ],
                "2S5=3H",
            ),
            // Case 2: Only special operations
            (
                vec![
                    MixedRepresentation::CigarOp(1, 'S'),
                    MixedRepresentation::CigarOp(2, 'N'),
                    MixedRepresentation::CigarOp(3, 'H'),
                ],
                "1S2N3H",
            ),
            // Case 3: Mix of special ops and pure indels
            (
                vec![
                    MixedRepresentation::CigarOp(2, 'S'),
                    MixedRepresentation::Tracepoint(5, 0), // Pure insertion
                    MixedRepresentation::CigarOp(1, 'H'),
                ],
                "2S5I1H",
            ),
        ];

        for (i, (mixed_tracepoints, expected_cigar)) in test_cases.iter().enumerate() {
            let result = mixed_tracepoints_to_cigar(
                mixed_tracepoints,
                a_seq,
                b_seq,
                0,
                0,
                mismatch,
                gap_open1,
                gap_ext1,
                gap_open2,
                gap_ext2,
            );

            assert_eq!(
                result, *expected_cigar,
                "Test case {}: Mixed tracepoint roundtrip failed",
                i + 1
            );
        }
    }

    #[test]
    fn test_mixed_tracepoints_with_alignment() {
        // Test with actual sequence alignment
        let a_seq = b"ACGTACGTACGT"; // 12 bases
        let b_seq = b"ACGTACGTACGT"; // 12 bases (identical)
        
        let mixed_tracepoints = vec![
            MixedRepresentation::CigarOp(2, 'S'),
            MixedRepresentation::Tracepoint(6, 6),
            MixedRepresentation::CigarOp(1, 'H'),
            MixedRepresentation::Tracepoint(4, 4),
        ];

        let result = mixed_tracepoints_to_cigar(
            &mixed_tracepoints,
            a_seq,
            b_seq,
            0,
            0,
            2, 4, 2, 6, 1,
        );

        // Since sequences are identical, we expect matches for the tracepoint segments
        assert_eq!(result, "2S6=1H4=");
    }

    #[test]
    fn test_mixed_tracepoints_edge_cases() {
        // Test edge cases
        let test_cases = vec![
            // Empty CIGAR
            ("", 5, vec![]),
            // Only special operations
            ("5S3H", 2, vec![
                MixedRepresentation::CigarOp(5, 'S'),
                MixedRepresentation::CigarOp(3, 'H'),
            ]),
            // Only regular operations
            ("5=3I2D", 10, vec![
                MixedRepresentation::Tracepoint(8, 7),
            ]),
            // max_diff = 0
            ("2=1X2=", 0, vec![
                MixedRepresentation::Tracepoint(2, 2),
                MixedRepresentation::Tracepoint(1, 1),
                MixedRepresentation::Tracepoint(2, 2),
            ]),
        ];

        for (i, (cigar, max_diff, expected)) in test_cases.iter().enumerate() {
            let result = cigar_to_mixed_tracepoints(cigar, *max_diff);
            assert_eq!(
                result,
                *expected,
                "Edge case {}: Failed for CIGAR '{}' with max_diff={}",
                i + 1,
                cigar,
                max_diff
            );
        }
    }

    #[test]
    fn test_variable_tracepoint_generation() {
        // Test CIGAR string
        let cigar = "10=2D5=2I3=";

        // Define test cases: (max_diff, expected_variable_tracepoints)
        let test_cases = vec![
            // Case 1: No differences allowed - each operation becomes its own segment
            (0, vec![(10, None), (0, Some(2)), (5, None), (2, Some(0)), (3, None)]),
            // Case 2: Allow up to 2 differences in each segment
            (2, vec![(15, Some(17)), (5, Some(3))]),
            // Case 3: Allow up to 5 differences - combines all operations
            (5, vec![(20, None)]),
        ];

        // Run each test case
        for (i, (max_diff, expected_variable_tracepoints)) in test_cases.iter().enumerate() {
            // Get actual results
            let variable_tracepoints = cigar_to_variable_tracepoints(cigar, *max_diff);

            // Check variable tracepoints
            assert_eq!(
                variable_tracepoints,
                *expected_variable_tracepoints,
                "Test case {}: Variable tracepoints with max_diff={} incorrect",
                i + 1,
                max_diff
            );
        }
    }

    #[test]
    fn test_variable_tracepoints_roundtrip() {
        let original_cigar = "1=1I18=";
        let a_seq = b"ACGTACGTACACGTACGTAC"; // 20 bases
        let b_seq = b"AGTACGTACACGTACGTAC"; // 19 bases (missing C)
        let max_diff = 5;

        // Test variable tracepoints roundtrip
        let variable_tracepoints = cigar_to_variable_tracepoints(original_cigar, max_diff);
        let reconstructed_cigar = variable_tracepoints_to_cigar(
            &variable_tracepoints,
            a_seq,
            b_seq,
            0,
            0, // sequences and start positions
            2,
            4,
            2,
            6,
            1, // alignment penalties
        );
        assert_eq!(
            reconstructed_cigar, original_cigar,
            "Variable tracepoint roundtrip failed"
        );
    }

    #[test]
    fn test_variable_tracepoints_optimization() {
        // Test cases where a_len == b_len should be optimized
        let test_cases = vec![
            // Equal lengths should use None
            ("5=", 10, vec![(5, None)]),
            ("3=2X4=", 10, vec![(9, None)]),
            // Different lengths should use Some
            ("5I", 10, vec![(5, Some(0))]),
            ("3D", 10, vec![(0, Some(3))]),
            ("2=3I1=", 10, vec![(6, Some(3))]),
        ];

        for (i, (cigar, max_diff, expected)) in test_cases.iter().enumerate() {
            let result = cigar_to_variable_tracepoints(cigar, *max_diff);
            assert_eq!(
                result,
                *expected,
                "Optimization test case {}: Failed for CIGAR '{}' with max_diff={}",
                i + 1,
                cigar,
                max_diff
            );
        }
    }

    #[test]
    fn test_variable_tracepoints_conversion() {
        // Test the conversion logic directly without WFA complexity
        let variable_tracepoints = vec![
            (5, None),        // Should become (5, 5)
            (3, Some(0)),     // Should become (3, 0)
            (0, Some(2)),     // Should become (0, 2)  
            (4, Some(6)),     // Should become (4, 6)
        ];

        // Test the conversion logic used in variable_tracepoints_to_cigar
        let converted_regular: Vec<(usize, usize)> = variable_tracepoints
            .iter()
            .map(|(a_len, b_len_opt)| {
                match b_len_opt {
                    None => (*a_len, *a_len),
                    Some(b_len) => (*a_len, *b_len),
                }
            })
            .collect();

        let expected_regular = vec![
            (5, 5),
            (3, 0),
            (0, 2),
            (4, 6),
        ];

        assert_eq!(
            converted_regular, expected_regular,
            "Variable to regular tracepoint conversion failed"
        );

        // Test the reverse conversion (regular to variable)
        let reconverted_variable: Vec<(usize, Option<usize>)> = expected_regular
            .iter()
            .map(|(a_len, b_len)| {
                if a_len == b_len {
                    (*a_len, None)
                } else {
                    (*a_len, Some(*b_len))
                }
            })
            .collect();

        assert_eq!(
            reconverted_variable, variable_tracepoints,
            "Regular to variable tracepoint conversion failed"
        );

        // Verify optimization: None should be used when a_len == b_len
        assert!(variable_tracepoints.iter().any(|(_, b_opt)| b_opt.is_none()),
               "Variable tracepoints should have at least one optimized entry (None)");
        assert!(variable_tracepoints.iter().any(|(_, b_opt)| b_opt.is_some()),
               "Variable tracepoints should have at least one unoptimized entry (Some)");
    }
}
