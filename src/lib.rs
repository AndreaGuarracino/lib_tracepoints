use lib_wfa2::affine_wavefront::{AffineWavefronts, AlignmentStatus};
use std::cmp::min;

/// Represents a CIGAR segment that can be either aligned or preserved as-is
#[derive(Debug, PartialEq)]
pub enum MixedRepresentation {
    /// Alignment segment represented by tracepoints
    Tracepoint(usize, usize),
    /// Special CIGAR operation that should be preserved intact
    CigarOp(usize, char),
}

/// Output type for tracepoint processing
enum TracepointOutput {
    /// Basic tracepoints as (a_len, b_len) pairs
    Basic(Vec<(usize, usize)>),
    /// Mixed representation with special CIGAR operations preserved
    Mixed(Vec<MixedRepresentation>),
}

/// Core function to process CIGAR into tracepoints with unified logic
/// 
/// Handles both basic and mixed tracepoint generation to eliminate code duplication.
/// Special operations (H, N, P, S) are preserved when preserve_special is true.
fn process_cigar_to_tracepoints(cigar: &str, max_diff: usize, preserve_special: bool) -> TracepointOutput {
    let ops = cigar_str_to_cigar_ops(cigar);
    let mut basic_tracepoints = Vec::new();
    let mut mixed_tracepoints = Vec::new();
    let mut cur_a_len = 0;
    let mut cur_b_len = 0;
    let mut cur_diff = 0;

    for (mut len, op) in ops {
        match op {
            'H' | 'N' | 'P' | 'S' if preserve_special => {
                // Preserve special operations as-is in mixed mode
                if cur_a_len > 0 || cur_b_len > 0 {
                    mixed_tracepoints.push(MixedRepresentation::Tracepoint(cur_a_len, cur_b_len));
                    cur_a_len = 0;
                    cur_b_len = 0;
                    cur_diff = 0;
                }
                mixed_tracepoints.push(MixedRepresentation::CigarOp(len, op));
            }
            'X' => {
                // Mismatches can be split across segments
                while len > 0 {
                    let remaining = max_diff.saturating_sub(cur_diff);
                    let step = min(len, remaining);

                    if step == 0 {
                        if cur_a_len > 0 || cur_b_len > 0 {
                            if preserve_special {
                                mixed_tracepoints.push(MixedRepresentation::Tracepoint(cur_a_len, cur_b_len));
                            } else {
                                basic_tracepoints.push((cur_a_len, cur_b_len));
                            }
                            cur_a_len = 0;
                            cur_b_len = 0;
                            cur_diff = 0;
                        }

                        if max_diff == 0 {
                            if preserve_special {
                                mixed_tracepoints.push(MixedRepresentation::Tracepoint(1, 1));
                            } else {
                                basic_tracepoints.push((1, 1));
                            }
                            len -= 1;
                        } else {
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
                            if preserve_special {
                                mixed_tracepoints.push(MixedRepresentation::Tracepoint(cur_a_len, cur_b_len));
                            } else {
                                basic_tracepoints.push((cur_a_len, cur_b_len));
                            }
                            cur_a_len = 0;
                            cur_b_len = 0;
                            cur_diff = 0;
                        }
                    }
                }
            }
            'I' | 'D' => {
                // Indels are kept together when possible
                if len > max_diff {
                    if cur_a_len > 0 || cur_b_len > 0 {
                        if preserve_special {
                            mixed_tracepoints.push(MixedRepresentation::Tracepoint(cur_a_len, cur_b_len));
                        } else {
                            basic_tracepoints.push((cur_a_len, cur_b_len));
                        }
                        cur_a_len = 0;
                        cur_b_len = 0;
                        cur_diff = 0;
                    }
                    let (a_add, b_add) = if op == 'I' { (len, 0) } else { (0, len) };
                    if preserve_special {
                        mixed_tracepoints.push(MixedRepresentation::Tracepoint(a_add, b_add));
                    } else {
                        basic_tracepoints.push((a_add, b_add));
                    }
                } else {
                    if cur_diff + len > max_diff {
                        if preserve_special {
                            mixed_tracepoints.push(MixedRepresentation::Tracepoint(cur_a_len, cur_b_len));
                        } else {
                            basic_tracepoints.push((cur_a_len, cur_b_len));
                        }
                        cur_a_len = 0;
                        cur_b_len = 0;
                        cur_diff = 0;
                    }
                    if op == 'I' {
                        cur_a_len += len;
                    } else {
                        cur_b_len += len;
                    }
                    cur_diff += len;
                }
            }
            '=' | 'M' => {
                // Matches don't count as differences
                cur_a_len += len;
                cur_b_len += len;
            }
            _ => panic!("Invalid CIGAR operation: {}", op),
        }
    }

    if cur_a_len > 0 || cur_b_len > 0 {
        if preserve_special {
            mixed_tracepoints.push(MixedRepresentation::Tracepoint(cur_a_len, cur_b_len));
        } else {
            basic_tracepoints.push((cur_a_len, cur_b_len));
        }
    }

    if preserve_special {
        TracepointOutput::Mixed(mixed_tracepoints)
    } else {
        TracepointOutput::Basic(basic_tracepoints)
    }
}

/// Convert CIGAR string into basic tracepoints
/// 
/// Segments CIGAR into tracepoints with at most max_diff differences per segment.
pub fn cigar_to_tracepoints(cigar: &str, max_diff: usize) -> Vec<(usize, usize)> {
    match process_cigar_to_tracepoints(cigar, max_diff, false) {
        TracepointOutput::Basic(tracepoints) => tracepoints,
        _ => unreachable!(),
    }
}

/// Convert CIGAR string into mixed representation tracepoints
/// 
/// Like cigar_to_tracepoints but preserves special operations (H, N, P, S).
pub fn cigar_to_mixed_tracepoints(cigar: &str, max_diff: usize) -> Vec<MixedRepresentation> {
    match process_cigar_to_tracepoints(cigar, max_diff, true) {
        TracepointOutput::Mixed(tracepoints) => tracepoints,
        _ => unreachable!(),
    }
}

/// Convert CIGAR string into variable tracepoints with length optimization
/// 
/// Uses (length, None) when a_len == b_len, otherwise (a_len, Some(b_len)).
pub fn cigar_to_variable_tracepoints(cigar: &str, max_diff: usize) -> Vec<(usize, Option<usize>)> {
    cigar_to_tracepoints(cigar, max_diff)
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

/// Reconstruct CIGAR string from tracepoint segments using WFA alignment
/// 
/// penalties: (mismatch, gap_open1, gap_ext1, gap_open2, gap_ext2)
fn reconstruct_cigar_from_segments(
    segments: &[(usize, usize)],
    a_seq: &[u8],
    b_seq: &[u8],
    a_start: usize,
    b_start: usize,
    penalties: (i32, i32, i32, i32, i32),
) -> String {
    let (mismatch, gap_open1, gap_ext1, gap_open2, gap_ext2) = penalties;
    let mut aligner = AffineWavefronts::with_penalties_affine2p(
        0, mismatch, gap_open1, gap_ext1, gap_open2, gap_ext2,
    );

    let mut cigar_ops = Vec::new();
    let mut current_a = a_start;
    let mut current_b = b_start;

    for &(a_len, b_len) in segments {
        if a_len > 0 && b_len == 0 {
            // Pure insertion
            cigar_ops.push((a_len, 'I'));
            current_a += a_len;
        } else if b_len > 0 && a_len == 0 {
            // Pure deletion
            cigar_ops.push((b_len, 'D'));
            current_b += b_len;
        } else {
            // Mixed segment - realign with WFA
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

/// Reconstruct CIGAR string from basic tracepoints
/// 
/// penalties: (mismatch, gap_open1, gap_ext1, gap_open2, gap_ext2)
pub fn tracepoints_to_cigar(
    tracepoints: &[(usize, usize)],
    a_seq: &[u8],
    b_seq: &[u8],
    a_start: usize,
    b_start: usize,
    penalties: (i32, i32, i32, i32, i32),
) -> String {
    reconstruct_cigar_from_segments(tracepoints, a_seq, b_seq, a_start, b_start, penalties)
}

/// Reconstruct CIGAR string from mixed representation tracepoints
/// 
/// Processes both alignment segments and preserved special operations.
/// penalties: (mismatch, gap_open1, gap_ext1, gap_open2, gap_ext2)
pub fn mixed_tracepoints_to_cigar(
    mixed_tracepoints: &[MixedRepresentation],
    a_seq: &[u8],
    b_seq: &[u8],
    a_start: usize,
    b_start: usize,
    penalties: (i32, i32, i32, i32, i32),
) -> String {
    let (mismatch, gap_open1, gap_ext1, gap_open2, gap_ext2) = penalties;
    let mut aligner = AffineWavefronts::with_penalties_affine2p(
        0, mismatch, gap_open1, gap_ext1, gap_open2, gap_ext2,
    );

    let mut cigar_ops = Vec::new();
    let mut current_a = a_start;
    let mut current_b = b_start;

    for item in mixed_tracepoints {
        match item {
            MixedRepresentation::CigarOp(len, op) => {
                // Add special operations directly
                cigar_ops.push((*len, *op));
            }
            MixedRepresentation::Tracepoint(a_len, b_len) => {
                // Process alignment segments
                if *a_len > 0 && *b_len == 0 {
                    cigar_ops.push((*a_len, 'I'));
                    current_a += a_len;
                } else if *b_len > 0 && *a_len == 0 {
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

/// Reconstruct CIGAR string from variable tracepoints
/// 
/// Converts variable format back to regular tracepoints, then reconstructs CIGAR.
/// penalties: (mismatch, gap_open1, gap_ext1, gap_open2, gap_ext2)
pub fn variable_tracepoints_to_cigar(
    variable_tracepoints: &[(usize, Option<usize>)],
    a_seq: &[u8],
    b_seq: &[u8],
    a_start: usize,
    b_start: usize,
    penalties: (i32, i32, i32, i32, i32),
) -> String {
    // Convert variable format: (len, None) -> (len, len), (a_len, Some(b_len)) -> (a_len, b_len)
    let regular_tracepoints: Vec<(usize, usize)> = variable_tracepoints
        .iter()
        .map(|(a_len, b_len_opt)| match b_len_opt {
            None => (*a_len, *a_len),
            Some(b_len) => (*a_len, *b_len),
        })
        .collect();
    
    reconstruct_cigar_from_segments(&regular_tracepoints, a_seq, b_seq, a_start, b_start, penalties)
}

// Helper functions

/// Merge consecutive CIGAR operations of the same type in-place
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

/// Parse CIGAR string into (length, operation) pairs
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

/// Convert WFA byte array to (length, operation) pairs
/// Treats 'M' (77) as '=' for consistency
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

    result.push((count, current_op));
    result
}

/// Format CIGAR operations as standard CIGAR string
pub fn cigar_ops_to_cigar_string(ops: &[(usize, char)]) -> String {
    ops.iter()
        .map(|(len, op)| format!("{}{}", len, op))
        .collect::<Vec<_>>()
        .join("")
}

/// Align two sequence segments using WFA algorithm
/// 
/// Note: aligns b vs a to get correct I/D orientation in CIGAR
pub fn align_sequences_wfa(
    a: &[u8],
    b: &[u8],
    aligner: &mut AffineWavefronts,
) -> Vec<(usize, char)> {
    let status = aligner.align(b, a);

    match status {
        AlignmentStatus::Completed => cigar_u8_to_cigar_ops(aligner.cigar()),
        s => panic!("Alignment failed with status: {:?}", s),
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
            0,
            (2, 4, 2, 6, 1),
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
                (mismatch, gap_open1, gap_ext1, gap_open2, gap_ext2),
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
            (2, 4, 2, 6, 1),
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
            0,
            (2, 4, 2, 6, 1),
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
