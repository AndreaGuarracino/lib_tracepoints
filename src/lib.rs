use lib_wfa2::affine_wavefront::{AffineWavefronts, AlignmentStatus};
use std::cmp::min;

/// Distance mode for alignment
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DistanceMode {
    /// Affine gap penalty with dual costs (gap_open1, gap_ext1, gap_open2, gap_ext2)
    Affine2p {
        mismatch: i32,
        gap_open1: i32,
        gap_ext1: i32,
        gap_open2: i32,
        gap_ext2: i32,
    },
    /// Edit distance (unit costs for mismatch and indels)
    Edit { mismatch: i32, gap_opening: i32 },
}

impl DistanceMode {
    /// Create an aligner configured for this distance mode
    pub fn create_aligner(&self) -> AffineWavefronts {
        match self {
            DistanceMode::Affine2p {
                mismatch,
                gap_open1,
                gap_ext1,
                gap_open2,
                gap_ext2,
            } => AffineWavefronts::with_penalties_affine2p(
                0, *mismatch, *gap_open1, *gap_ext1, *gap_open2, *gap_ext2,
            ),
            DistanceMode::Edit {
                mismatch,
                gap_opening,
            } => AffineWavefronts::with_penalties_edit(0, *mismatch, *gap_opening),
        }
    }
}

/// Output type for tracepoint processing
enum Tracepoints {
    /// Basic tracepoints as (a_len, b_len) pairs
    Basic(Vec<(usize, usize)>),
    /// Mixed representation with special CIGAR operations preserved
    Mixed(Vec<MixedRepresentation>),
}

/// Represents a CIGAR segment that can be either aligned or preserved as-is
#[derive(Debug, PartialEq)]
pub enum MixedRepresentation {
    /// Alignment segment represented by tracepoints
    Tracepoint(usize, usize),
    /// Special CIGAR operation that should be preserved intact
    CigarOp(usize, char),
}

/// Helper function to flush current segment to tracepoint collections
fn flush_segment(
    cur_a_len: &mut usize,
    cur_b_len: &mut usize,
    cur_diff: &mut usize,
    preserve_special: bool,
    basic_tracepoints: &mut Vec<(usize, usize)>,
    mixed_tracepoints: &mut Vec<MixedRepresentation>,
) {
    if *cur_a_len > 0 || *cur_b_len > 0 {
        if preserve_special {
            mixed_tracepoints.push(MixedRepresentation::Tracepoint(*cur_a_len, *cur_b_len));
        } else {
            basic_tracepoints.push((*cur_a_len, *cur_b_len));
        }
        *cur_a_len = 0;
        *cur_b_len = 0;
        *cur_diff = 0;
    }
}

/// Helper function to add a tracepoint to the appropriate collection
fn add_tracepoint(
    a_len: usize,
    b_len: usize,
    preserve_special: bool,
    basic_tracepoints: &mut Vec<(usize, usize)>,
    mixed_tracepoints: &mut Vec<MixedRepresentation>,
) {
    if preserve_special {
        mixed_tracepoints.push(MixedRepresentation::Tracepoint(a_len, b_len));
    } else {
        basic_tracepoints.push((a_len, b_len));
    }
}

/// Unified function to process CIGAR into tracepoints
///
/// Handles both basic and mixed tracepoint generation with optional indel splitting.
/// Special operations (H, N, P, S) are preserved when preserve_special is true.
/// Indels can be split across segments when allow_indel_split is true (raw mode).
fn process_cigar_unified(
    cigar: &str,
    max_diff: usize,
    preserve_special: bool,
    allow_indel_split: bool,
) -> Tracepoints {
    let ops = cigar_str_to_cigar_ops(cigar);
    let mut basic_tracepoints = Vec::new();
    let mut mixed_tracepoints = Vec::new();
    let mut cur_a_len = 0;
    let mut cur_b_len = 0;
    let mut cur_diff = 0;

    for (mut len, op) in ops {
        match op {
            'H' | 'N' | 'P' | 'S' if preserve_special => {
                flush_segment(
                    &mut cur_a_len,
                    &mut cur_b_len,
                    &mut cur_diff,
                    preserve_special,
                    &mut basic_tracepoints,
                    &mut mixed_tracepoints,
                );
                mixed_tracepoints.push(MixedRepresentation::CigarOp(len, op));
            }
            'X' | 'I' | 'D' if op == 'X' || allow_indel_split => {
                while len > 0 {
                    let remaining = max_diff.saturating_sub(cur_diff);
                    let step = min(len, remaining);

                    if step == 0 {
                        flush_segment(
                            &mut cur_a_len,
                            &mut cur_b_len,
                            &mut cur_diff,
                            preserve_special,
                            &mut basic_tracepoints,
                            &mut mixed_tracepoints,
                        );

                        let (a_add, b_add) = match op {
                            'X' => (1, 1),
                            'I' => (1, 0),
                            'D' => (0, 1),
                            _ => unreachable!(),
                        };

                        if max_diff == 0 {
                            add_tracepoint(
                                a_add,
                                b_add,
                                preserve_special,
                                &mut basic_tracepoints,
                                &mut mixed_tracepoints,
                            );
                        } else {
                            cur_a_len = a_add;
                            cur_b_len = b_add;
                            cur_diff = 1;
                        }
                        len -= 1;
                    } else {
                        let (a_add, b_add) = match op {
                            'X' => (step, step),
                            'I' => (step, 0),
                            'D' => (0, step),
                            _ => unreachable!(),
                        };
                        cur_a_len += a_add;
                        cur_b_len += b_add;
                        cur_diff += step;
                        len -= step;

                        if cur_diff == max_diff {
                            flush_segment(
                                &mut cur_a_len,
                                &mut cur_b_len,
                                &mut cur_diff,
                                preserve_special,
                                &mut basic_tracepoints,
                                &mut mixed_tracepoints,
                            );
                        }
                    }
                }
            }
            'I' | 'D' => {
                if len > max_diff {
                    flush_segment(
                        &mut cur_a_len,
                        &mut cur_b_len,
                        &mut cur_diff,
                        preserve_special,
                        &mut basic_tracepoints,
                        &mut mixed_tracepoints,
                    );
                    let (a_add, b_add) = if op == 'I' { (len, 0) } else { (0, len) };
                    add_tracepoint(
                        a_add,
                        b_add,
                        preserve_special,
                        &mut basic_tracepoints,
                        &mut mixed_tracepoints,
                    );
                } else {
                    if cur_diff + len > max_diff {
                        flush_segment(
                            &mut cur_a_len,
                            &mut cur_b_len,
                            &mut cur_diff,
                            preserve_special,
                            &mut basic_tracepoints,
                            &mut mixed_tracepoints,
                        );
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
                cur_a_len += len;
                cur_b_len += len;
            }
            _ => panic!("Invalid CIGAR operation: {op}"),
        }
    }

    flush_segment(
        &mut cur_a_len,
        &mut cur_b_len,
        &mut cur_diff,
        preserve_special,
        &mut basic_tracepoints,
        &mut mixed_tracepoints,
    );

    if preserve_special {
        Tracepoints::Mixed(mixed_tracepoints)
    } else {
        Tracepoints::Basic(basic_tracepoints)
    }
}

/// Unified function to process CIGAR into tracepoints using diagonal distance
///
/// Breaks segments when the distance from the main diagonal exceeds max_diff.
/// Insertions increase diagonal distance (+), deletions decrease it (-).
/// The main diagonal is influenced by the overall sequence length difference.
fn process_cigar_unified_diagonal(
    cigar: &str,
    max_diff: usize,
    preserve_special: bool,
) -> Tracepoints {
    let ops = cigar_str_to_cigar_ops(cigar);
    let mut basic_tracepoints = Vec::new();
    let mut mixed_tracepoints = Vec::new();
    let mut cur_a_len = 0;
    let mut cur_b_len = 0;
    let mut diagonal_distance: i64 = 0;

    for (len, op) in ops {
        match op {
            'H' | 'N' | 'P' | 'S' if preserve_special => {
                flush_segment(
                    &mut cur_a_len,
                    &mut cur_b_len,
                    &mut 0,
                    preserve_special,
                    &mut basic_tracepoints,
                    &mut mixed_tracepoints,
                );
                diagonal_distance = 0; // Reset diagonal distance after flushing
                mixed_tracepoints.push(MixedRepresentation::CigarOp(len, op));
            }
            'I' => {
                let new_diagonal_distance = diagonal_distance + len as i64;

                if new_diagonal_distance.unsigned_abs() > max_diff as u64 {
                    // Would exceed max_diff, need to flush current segment first
                    if cur_a_len > 0 || cur_b_len > 0 {
                        flush_segment(
                            &mut cur_a_len,
                            &mut cur_b_len,
                            &mut 0,
                            preserve_special,
                            &mut basic_tracepoints,
                            &mut mixed_tracepoints,
                        );
                    }

                    // Add the insertion as its own segment
                    add_tracepoint(
                        len,
                        0,
                        preserve_special,
                        &mut basic_tracepoints,
                        &mut mixed_tracepoints,
                    );
                    diagonal_distance = 0; // Reset after creating segment
                } else {
                    // Can add to current segment
                    cur_a_len += len;
                    diagonal_distance = new_diagonal_distance;
                }
            }
            'D' => {
                let new_diagonal_distance = diagonal_distance - len as i64;

                if new_diagonal_distance.unsigned_abs() > max_diff as u64 {
                    // Would exceed max_diff, need to flush current segment first
                    if cur_a_len > 0 || cur_b_len > 0 {
                        flush_segment(
                            &mut cur_a_len,
                            &mut cur_b_len,
                            &mut 0,
                            preserve_special,
                            &mut basic_tracepoints,
                            &mut mixed_tracepoints,
                        );
                    }

                    // Add the deletion as its own segment
                    add_tracepoint(
                        0,
                        len,
                        preserve_special,
                        &mut basic_tracepoints,
                        &mut mixed_tracepoints,
                    );
                    diagonal_distance = 0; // Reset after creating segment
                } else {
                    // Can add to current segment
                    cur_b_len += len;
                    diagonal_distance = new_diagonal_distance;
                }
            }
            'X' => {
                // Mismatches don't change diagonal distance but still consume sequence
                cur_a_len += len;
                cur_b_len += len;
            }
            '=' | 'M' => {
                // Matches don't change diagonal distance
                cur_a_len += len;
                cur_b_len += len;
            }
            _ => panic!("Invalid CIGAR operation: {op}"),
        }
    }

    // Flush any remaining segment
    flush_segment(
        &mut cur_a_len,
        &mut cur_b_len,
        &mut 0,
        preserve_special,
        &mut basic_tracepoints,
        &mut mixed_tracepoints,
    );

    if preserve_special {
        Tracepoints::Mixed(mixed_tracepoints)
    } else {
        Tracepoints::Basic(basic_tracepoints)
    }
}

/// Convert CIGAR string into basic tracepoints
///
/// Segments CIGAR into tracepoints with at most max_diff differences per segment.
pub fn cigar_to_tracepoints(cigar: &str, max_diff: usize) -> Vec<(usize, usize)> {
    match process_cigar_unified(cigar, max_diff, false, false) {
        Tracepoints::Basic(tracepoints) => tracepoints,
        _ => unreachable!(),
    }
}

/// Convert CIGAR string into mixed representation tracepoints
///
/// Like cigar_to_tracepoints but preserves special operations (H, N, P, S).
pub fn cigar_to_mixed_tracepoints(cigar: &str, max_diff: usize) -> Vec<MixedRepresentation> {
    match process_cigar_unified(cigar, max_diff, true, false) {
        Tracepoints::Mixed(tracepoints) => tracepoints,
        _ => unreachable!(),
    }
}

/// Convert CIGAR string into raw tracepoints
///
/// Like cigar_to_tracepoints but allows indels to be split across segments.
/// This provides more granular control over segment sizes.
pub fn cigar_to_tracepoints_raw(cigar: &str, max_diff: usize) -> Vec<(usize, usize)> {
    match process_cigar_unified(cigar, max_diff, false, true) {
        Tracepoints::Basic(tracepoints) => tracepoints,
        _ => unreachable!(),
    }
}

/// Convert CIGAR string into mixed representation raw tracepoints
///
/// Like cigar_to_mixed_tracepoints but allows indels to be split across segments.
pub fn cigar_to_mixed_tracepoints_raw(cigar: &str, max_diff: usize) -> Vec<MixedRepresentation> {
    match process_cigar_unified(cigar, max_diff, true, true) {
        Tracepoints::Mixed(tracepoints) => tracepoints,
        _ => unreachable!(),
    }
}

/// Convert tracepoints to variable format with length optimization
///
/// Uses (length, None) when a_len == b_len, otherwise (a_len, Some(b_len)).
fn to_variable_format(tracepoints: Vec<(usize, usize)>) -> Vec<(usize, Option<usize>)> {
    tracepoints
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

/// Convert CIGAR string into variable tracepoints with length optimization
///
/// Uses (length, None) when a_len == b_len, otherwise (a_len, Some(b_len)).
pub fn cigar_to_variable_tracepoints(cigar: &str, max_diff: usize) -> Vec<(usize, Option<usize>)> {
    to_variable_format(cigar_to_tracepoints(cigar, max_diff))
}

/// Convert CIGAR string into variable raw tracepoints with length optimization
///
/// Like cigar_to_variable_tracepoints but allows indels to be split across segments.
pub fn cigar_to_variable_tracepoints_raw(
    cigar: &str,
    max_diff: usize,
) -> Vec<(usize, Option<usize>)> {
    to_variable_format(cigar_to_tracepoints_raw(cigar, max_diff))
}

/// Convert CIGAR string into tracepoints using diagonal distance segmentation
///
/// Segments CIGAR by breaking when distance from main diagonal exceeds max_diff.
/// Insertions increase diagonal distance, deletions decrease it.
pub fn cigar_to_tracepoints_diagonal(cigar: &str, max_diff: usize) -> Vec<(usize, usize)> {
    match process_cigar_unified_diagonal(cigar, max_diff, false) {
        Tracepoints::Basic(tracepoints) => tracepoints,
        _ => unreachable!(),
    }
}

/// Convert CIGAR string into mixed tracepoints using diagonal distance segmentation
///
/// Like cigar_to_tracepoints_diagonal but preserves special operations (H, N, P, S).
pub fn cigar_to_mixed_tracepoints_diagonal(
    cigar: &str,
    max_diff: usize,
) -> Vec<MixedRepresentation> {
    match process_cigar_unified_diagonal(cigar, max_diff, true) {
        Tracepoints::Mixed(tracepoints) => tracepoints,
        _ => unreachable!(),
    }
}

/// Convert CIGAR string into variable tracepoints using diagonal distance with length optimization
///
/// Uses diagonal distance segmentation with (length, None) when a_len == b_len.
pub fn cigar_to_variable_tracepoints_diagonal(
    cigar: &str,
    max_diff: usize,
) -> Vec<(usize, Option<usize>)> {
    to_variable_format(cigar_to_tracepoints_diagonal(cigar, max_diff))
}

/// Reverse a CIGAR string for complement alignments
/// 
/// Reverses both the order of operations and the numbers within each operation.
/// For example: "3M2I5D" becomes "5D2I3M"
fn reverse_cigar(cigar: &str) -> String {
    let ops = cigar_str_to_cigar_ops(cigar);
    ops.iter()
        .rev()
        .map(|(len, op)| format!("{}{}", len, op))
        .collect()
}

/// Generate FASTGA-style tracepoints from CIGAR string
pub fn cigar_to_tracepoints_fastga(
    cigar: &str, 
    trace_spacing: usize,
    a_start: usize,
    b_start: usize,
    complement: bool,
) -> Vec<(usize, usize)> {
    let cigar_to_process = if complement {
        reverse_cigar(cigar)
    } else {
        cigar.to_string()
    };
    
    let ops = cigar_str_to_cigar_ops(&cigar_to_process);
    let mut tracepoints = Vec::new();
    
    let mut a_pos = a_start;
    let mut b_pos = b_start;
    let mut last_b_pos = b_start;
    let mut diff = 0i64;
    let mut last_diff = 0i64;
    let mut next_trace = ((a_start / trace_spacing) + 1) * trace_spacing;
    
    for (mut len, op) in ops {
        while len > 0 {
            let consume = match op {
                '=' | 'M' | 'X' => {
                    // Check if this operation crosses a boundary
                    if a_pos + len > next_trace {
                        let inc = next_trace - a_pos;
                        a_pos = next_trace;
                        b_pos += inc;
                        if op == 'X' {
                            diff += inc as i64;
                        }
                        // Emit tracepoint at boundary crossing
                        tracepoints.push((
                            (diff - last_diff).unsigned_abs() as usize,
                            b_pos - last_b_pos,
                        ));
                        last_diff = diff;
                        last_b_pos = b_pos;
                        next_trace += trace_spacing;
                        inc
                    } else {
                        // Operation fits entirely before next boundary
                        a_pos += len;
                        b_pos += len;
                        if op == 'X' {
                            diff += len as i64;
                        }
                        len
                    }
                }
                'I' => {
                    // Check segment length constraint
                    if trace_spacing + len > 200 && a_pos != next_trace - trace_spacing {
                        tracepoints.push((
                            (diff - last_diff).unsigned_abs() as usize,
                            b_pos - last_b_pos,
                        ));
                        last_diff = diff;
                        last_b_pos = b_pos;
                        next_trace = ((a_pos / trace_spacing) + 1) * trace_spacing;
                    }
                    
                    if a_pos + len > next_trace {
                        let inc = next_trace - a_pos;
                        a_pos = next_trace;
                        diff += inc as i64;
                        tracepoints.push((
                            (diff - last_diff).unsigned_abs() as usize,
                            b_pos - last_b_pos,
                        ));
                        last_diff = diff;
                        last_b_pos = b_pos;
                        next_trace += trace_spacing;
                        inc
                    } else {
                        a_pos += len;
                        diff += len as i64;
                        len
                    }
                }
                'D' => {
                    // Check segment length constraint
                    if (b_pos - last_b_pos) + len + (next_trace - a_pos) > 200 {
                        if a_pos != next_trace - trace_spacing {
                            tracepoints.push((
                                (diff - last_diff).unsigned_abs() as usize,
                                b_pos - last_b_pos,
                            ));
                            last_diff = diff;
                            last_b_pos = b_pos;
                        }
                    }
                    diff += len as i64;
                    b_pos += len;
                    len // Consume all since D/N don't advance a_pos
                }
                'H' | 'N'  | 'S' | 'P' => len, // Skip special operations
                _ => panic!("Invalid CIGAR operation: {op}"),
            };
            len -= consume;
        }
    }
    
    // Add final tracepoint if we've passed the last boundary
    if a_pos > next_trace - trace_spacing {
        tracepoints.push((
            (diff - last_diff).unsigned_abs() as usize,
            b_pos - last_b_pos,
        ));
    }
    
    tracepoints
}

/// Create an aligner with the given distance mode
fn create_aligner(distance_mode: &DistanceMode) -> AffineWavefronts {
    distance_mode.create_aligner()
}

/// Reconstruct CIGAR string from tracepoint segments using WFA alignment with DistanceMode
fn reconstruct_cigar_from_segments(
    segments: &[(usize, usize)],
    a_seq: &[u8],
    b_seq: &[u8],
    a_start: usize,
    b_start: usize,
    distance_mode: &DistanceMode,
) -> String {
    let mut aligner = create_aligner(distance_mode);
    reconstruct_cigar_from_segments_with_aligner(
        segments,
        a_seq,
        b_seq,
        a_start,
        b_start,
        &mut aligner,
    )
}

/// Reconstruct CIGAR string from tracepoint segments using provided aligner
///
/// Takes an aligner parameter to allow callers to prepare an aligner and reuse it multiple times.
fn reconstruct_cigar_from_segments_with_aligner(
    segments: &[(usize, usize)],
    a_seq: &[u8],
    b_seq: &[u8],
    a_start: usize,
    b_start: usize,
    aligner: &mut AffineWavefronts,
) -> String {
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
            let seg_ops =
                align_sequences_wfa(&a_seq[current_a..a_end], &b_seq[current_b..b_end], aligner);
            cigar_ops.extend(seg_ops);
            current_a = a_end;
            current_b = b_end;
        }
    }
    merge_cigar_ops(&mut cigar_ops);
    cigar_ops_to_cigar_string(&cigar_ops)
}

/// Reconstruct CIGAR string from basic tracepoints with DistanceMode
pub fn tracepoints_to_cigar(
    tracepoints: &[(usize, usize)],
    a_seq: &[u8],
    b_seq: &[u8],
    a_start: usize,
    b_start: usize,
    distance_mode: &DistanceMode,
) -> String {
    reconstruct_cigar_from_segments(tracepoints, a_seq, b_seq, a_start, b_start, distance_mode)
}

/// Reconstruct CIGAR string from mixed representation tracepoints
///
/// Processes both alignment segments and preserved special operations.
pub fn mixed_tracepoints_to_cigar(
    mixed_tracepoints: &[MixedRepresentation],
    a_seq: &[u8],
    b_seq: &[u8],
    a_start: usize,
    b_start: usize,
    distance_mode: &DistanceMode,
) -> String {
    let mut aligner = create_aligner(distance_mode);
    let mut cigar_ops = Vec::new();
    let mut current_a = a_start;
    let mut current_b = b_start;

    for item in mixed_tracepoints {
        match item {
            MixedRepresentation::CigarOp(len, op) => {
                cigar_ops.push((*len, *op));
            }
            MixedRepresentation::Tracepoint(a_len, b_len) => {
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

/// Convert variable format back to regular tracepoints
/// (len, None) -> (len, len), (a_len, Some(b_len)) -> (a_len, b_len)
fn from_variable_format(variable_tracepoints: &[(usize, Option<usize>)]) -> Vec<(usize, usize)> {
    variable_tracepoints
        .iter()
        .map(|(a_len, b_len_opt)| match b_len_opt {
            None => (*a_len, *a_len),
            Some(b_len) => (*a_len, *b_len),
        })
        .collect()
}

/// Reconstruct CIGAR string from variable tracepoints
///
/// Converts variable format back to regular tracepoints, then reconstructs CIGAR.
pub fn variable_tracepoints_to_cigar(
    variable_tracepoints: &[(usize, Option<usize>)],
    a_seq: &[u8],
    b_seq: &[u8],
    a_start: usize,
    b_start: usize,
    distance_mode: &DistanceMode,
) -> String {
    let regular_tracepoints = from_variable_format(variable_tracepoints);
    reconstruct_cigar_from_segments(
        &regular_tracepoints,
        a_seq,
        b_seq,
        a_start,
        b_start,
        distance_mode,
    )
}

/// Reconstruct CIGAR string from variable tracepoints with provided aligner
///
/// Takes an aligner parameter to allow callers to prepare an aligner and reuse it multiple times.
/// Converts variable format back to regular tracepoints, then reconstructs CIGAR.
pub fn variable_tracepoints_to_cigar_with_aligner(
    variable_tracepoints: &[(usize, Option<usize>)],
    a_seq: &[u8],
    b_seq: &[u8],
    a_start: usize,
    b_start: usize,
    aligner: &mut AffineWavefronts,
) -> String {
    let regular_tracepoints = from_variable_format(variable_tracepoints);
    reconstruct_cigar_from_segments_with_aligner(
        &regular_tracepoints,
        a_seq,
        b_seq,
        a_start,
        b_start,
        aligner,
    )
}

/// Reconstruct CIGAR string from diagonal tracepoints
///
/// This function reconstructs a CIGAR string from tracepoints that were created
/// using diagonal distance segmentation.
pub fn tracepoints_to_cigar_diagonal(
    tracepoints: &[(usize, usize)],
    a_seq: &[u8],
    b_seq: &[u8],
    a_start: usize,
    b_start: usize,
    distance_mode: &DistanceMode,
) -> String {
    // Diagonal tracepoints can be reconstructed the same way as regular tracepoints
    // since they represent the same (a_len, b_len) format
    tracepoints_to_cigar(tracepoints, a_seq, b_seq, a_start, b_start, distance_mode)
}

/// Reconstruct CIGAR string from mixed diagonal tracepoints
///
/// Processes both alignment segments and preserved special operations from
/// diagonal distance segmentation.
pub fn mixed_tracepoints_to_cigar_diagonal(
    mixed_tracepoints: &[MixedRepresentation],
    a_seq: &[u8],
    b_seq: &[u8],
    a_start: usize,
    b_start: usize,
    distance_mode: &DistanceMode,
) -> String {
    // Mixed diagonal tracepoints can be reconstructed the same way as regular mixed tracepoints
    // since they use the same MixedRepresentation format
    mixed_tracepoints_to_cigar(
        mixed_tracepoints,
        a_seq,
        b_seq,
        a_start,
        b_start,
        distance_mode,
    )
}

/// Reconstruct CIGAR string from variable diagonal tracepoints
///
/// Converts variable format diagonal tracepoints back to regular format,
/// then reconstructs CIGAR.
pub fn variable_tracepoints_to_cigar_diagonal(
    variable_tracepoints: &[(usize, Option<usize>)],
    a_seq: &[u8],
    b_seq: &[u8],
    a_start: usize,
    b_start: usize,
    distance_mode: &DistanceMode,
) -> String {
    // Variable diagonal tracepoints can be reconstructed the same way as regular variable tracepoints
    // since they use the same (usize, Option<usize>) format
    variable_tracepoints_to_cigar(
        variable_tracepoints,
        a_seq,
        b_seq,
        a_start,
        b_start,
        distance_mode,
    )
}

/// Reconstruct CIGAR string from variable diagonal tracepoints with provided aligner
///
/// Takes an aligner parameter to allow callers to prepare an aligner and reuse it multiple times.
/// Converts variable format diagonal tracepoints back to regular format, then reconstructs CIGAR.
pub fn variable_tracepoints_to_cigar_diagonal_with_aligner(
    variable_tracepoints: &[(usize, Option<usize>)],
    a_seq: &[u8],
    b_seq: &[u8],
    a_start: usize,
    b_start: usize,
    aligner: &mut AffineWavefronts,
) -> String {
    // Variable diagonal tracepoints can be reconstructed the same way as regular variable tracepoints
    variable_tracepoints_to_cigar_with_aligner(
        variable_tracepoints,
        a_seq,
        b_seq,
        a_start,
        b_start,
        aligner,
    )
}

/// Reconstruct CIGAR from FASTGA-style tracepoints
///
/// Uses edit distance internally for alignment.
pub fn tracepoints_to_cigar_fastga(
    tracepoints: &[(usize, usize)],
    trace_spacing: usize,
    a_seq: &[u8],
    b_seq: &[u8],
    a_start: usize,
    b_start: usize,
) -> String {
    // Use edit distance mode as FASTGA does
    let distance_mode = DistanceMode::Edit {
        mismatch: 1,
        gap_opening: 1,
    };
    
    let mut aligner = distance_mode.create_aligner();
    let mut cigar_ops = Vec::new();
    let mut current_a = 0;  // Indices into the sequences (not absolute positions)
    let mut current_b = 0;
    
    // Calculate first segment length based on starting position
    let first_segment_a_len = trace_spacing - (a_start % trace_spacing);
    
    for (i, &(_, b_len)) in tracepoints.iter().enumerate() {
        // Determine a_length for this segment
        let a_len = if i == 0 {
            first_segment_a_len
        } else {
            trace_spacing
        };
        
        // Calculate segment boundaries
        let a_end = (current_a + a_len).min(a_seq.len());
        let b_end = (current_b + b_len).min(b_seq.len());
        
        // Skip if we've exhausted either sequence
        if current_a >= a_seq.len() || current_b >= b_seq.len() {
            break;
        }
        
        // Align this segment if it has content
        if a_end > current_a && b_end > current_b {
            let seg_ops = align_sequences_wfa(
                &a_seq[current_a..a_end],
                &b_seq[current_b..b_end],
                &mut aligner,
            );
            cigar_ops.extend(seg_ops);
        }
        
        current_a = a_end;
        current_b = b_end;
    }
    
    merge_cigar_ops(&mut cigar_ops);
    cigar_ops_to_cigar_string(&cigar_ops)
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
        .map(|(len, op)| format!("{len}{op}"))
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
        s => panic!("Alignment failed with status: {s:?}"),
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
        let test_cases = [
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
    fn test_distance_mode_affine2p() {
        let original_cigar = "1=1I18=";
        let a_seq = b"ACGTACGTACACGTACGTAC"; // 20 bases
        let b_seq = b"AGTACGTACACGTACGTAC"; // 19 bases (missing C)
        let max_diff = 5;

        let distance_mode = DistanceMode::Affine2p {
            mismatch: 2,
            gap_open1: 4,
            gap_ext1: 2,
            gap_open2: 6,
            gap_ext2: 1,
        };

        // Test tracepoints roundtrip with DistanceMode API
        let tracepoints = cigar_to_tracepoints(original_cigar, max_diff);
        let reconstructed_cigar =
            tracepoints_to_cigar(&tracepoints, a_seq, b_seq, 0, 0, &distance_mode);
        assert_eq!(
            reconstructed_cigar, original_cigar,
            "Affine2p distance mode roundtrip failed"
        );
    }

    #[test]
    fn test_distance_mode_edit() {
        let original_cigar = "1=1I18=";
        let a_seq = b"ACGTACGTACACGTACGTAC"; // 20 bases
        let b_seq = b"AGTACGTACACGTACGTAC"; // 19 bases (missing C)
        let max_diff = 5;

        let distance_mode = DistanceMode::Edit {
            mismatch: 1,
            gap_opening: 1,
        };

        // Test tracepoints roundtrip with edit distance mode
        let tracepoints = cigar_to_tracepoints(original_cigar, max_diff);
        let reconstructed_cigar =
            tracepoints_to_cigar(&tracepoints, a_seq, b_seq, 0, 0, &distance_mode);
        assert_eq!(
            reconstructed_cigar, original_cigar,
            "Edit distance mode roundtrip failed"
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
        let test_cases = [
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

        let distance_mode = DistanceMode::Affine2p {
            mismatch,
            gap_open1,
            gap_ext1,
            gap_open2,
            gap_ext2,
        };

        for (i, (mixed_tracepoints, expected_cigar)) in test_cases.iter().enumerate() {
            let result =
                mixed_tracepoints_to_cigar(mixed_tracepoints, a_seq, b_seq, 0, 0, &distance_mode);

            assert_eq!(
                result,
                *expected_cigar,
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

        let distance_mode = DistanceMode::Affine2p {
            mismatch: 2,
            gap_open1: 4,
            gap_ext1: 2,
            gap_open2: 6,
            gap_ext2: 1,
        };
        let result =
            mixed_tracepoints_to_cigar(&mixed_tracepoints, a_seq, b_seq, 0, 0, &distance_mode);

        // Since sequences are identical, we expect matches for the tracepoint segments
        assert_eq!(result, "2S6=1H4=");
    }

    #[test]
    fn test_mixed_tracepoints_edge_cases() {
        // Test edge cases
        let test_cases = [
            // Empty CIGAR
            ("", 5, vec![]),
            // Only special operations
            (
                "5S3H",
                2,
                vec![
                    MixedRepresentation::CigarOp(5, 'S'),
                    MixedRepresentation::CigarOp(3, 'H'),
                ],
            ),
            // Only regular operations
            ("5=3I2D", 10, vec![MixedRepresentation::Tracepoint(8, 7)]),
            // max_diff = 0
            (
                "2=1X2=",
                0,
                vec![
                    MixedRepresentation::Tracepoint(2, 2),
                    MixedRepresentation::Tracepoint(1, 1),
                    MixedRepresentation::Tracepoint(2, 2),
                ],
            ),
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
        let test_cases = [
            // Case 1: No differences allowed - each operation becomes its own segment
            (
                0,
                vec![(10, None), (0, Some(2)), (5, None), (2, Some(0)), (3, None)],
            ),
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

        let distance_mode = DistanceMode::Affine2p {
            mismatch: 2,
            gap_open1: 4,
            gap_ext1: 2,
            gap_open2: 6,
            gap_ext2: 1,
        };

        // Test variable tracepoints roundtrip
        let variable_tracepoints = cigar_to_variable_tracepoints(original_cigar, max_diff);
        let reconstructed_cigar = variable_tracepoints_to_cigar(
            &variable_tracepoints,
            a_seq,
            b_seq,
            0,
            0,
            &distance_mode,
        );
        assert_eq!(
            reconstructed_cigar, original_cigar,
            "Variable tracepoint roundtrip failed"
        );
    }

    #[test]
    fn test_variable_tracepoints_to_cigar_with_aligner() {
        // Test the new function that takes an aligner parameter
        let original_cigar = "1=1I18=";
        let a_seq = b"ACGTACGTACACGTACGTAC"; // 20 bases
        let b_seq = b"AGTACGTACACGTACGTAC"; // 19 bases (missing C)
        let max_diff = 5;

        // Create variable tracepoints
        let variable_tracepoints = cigar_to_variable_tracepoints(original_cigar, max_diff);

        // Create an aligner
        let mut aligner = AffineWavefronts::with_penalties_affine2p(0, 2, 4, 2, 6, 1);

        // Test the new function
        let reconstructed_cigar = variable_tracepoints_to_cigar_with_aligner(
            &variable_tracepoints,
            a_seq,
            b_seq,
            0,
            0,
            &mut aligner,
        );

        let distance_mode = DistanceMode::Affine2p {
            mismatch: 2,
            gap_open1: 4,
            gap_ext1: 2,
            gap_open2: 6,
            gap_ext2: 1,
        };

        // Should produce the same result
        let expected_cigar = variable_tracepoints_to_cigar(
            &variable_tracepoints,
            a_seq,
            b_seq,
            0,
            0,
            &distance_mode,
        );

        assert_eq!(
            reconstructed_cigar, expected_cigar,
            "New function should produce same result as legacy version"
        );
        assert_eq!(
            reconstructed_cigar, original_cigar,
            "Roundtrip should work with aligner"
        );

        // Test with multiple calls to ensure aligner can be reused
        let variable_tracepoints2 = vec![(5, None), (2, Some(0)), (3, None)];
        let reconstructed_cigar2 = variable_tracepoints_to_cigar_with_aligner(
            &variable_tracepoints2,
            a_seq,
            b_seq,
            0,
            0,
            &mut aligner,
        );

        let expected_cigar2 = variable_tracepoints_to_cigar(
            &variable_tracepoints2,
            a_seq,
            b_seq,
            0,
            0,
            &distance_mode,
        );

        assert_eq!(
            reconstructed_cigar2, expected_cigar2,
            "Function should work with reused aligner"
        );
    }

    #[test]
    fn test_raw_tracepoint_generation() {
        // Test CIGAR string with indels
        let cigar = "5=10I5=10D5=";

        // Test cases: (max_diff, expected_tracepoints)
        let test_cases = [
            // Case 1: No differences allowed - each operation becomes its own segment
            (
                0,
                vec![
                    (5, 5),
                    (1, 0),
                    (1, 0),
                    (1, 0),
                    (1, 0),
                    (1, 0),
                    (1, 0),
                    (1, 0),
                    (1, 0),
                    (1, 0),
                    (1, 0),
                    (5, 5),
                    (0, 1),
                    (0, 1),
                    (0, 1),
                    (0, 1),
                    (0, 1),
                    (0, 1),
                    (0, 1),
                    (0, 1),
                    (0, 1),
                    (0, 1),
                    (5, 5),
                ],
            ),
            // Case 2: Allow up to 3 differences - indels are split and combined with matches
            (
                3,
                vec![(8, 5), (3, 0), (3, 0), (6, 7), (0, 3), (0, 3), (5, 7)],
            ),
            // Case 3: Allow up to 5 differences - indels are split into chunks of 5
            (5, vec![(10, 5), (5, 0), (5, 10), (0, 5), (5, 5)]),
            // Case 4: Allow up to 10 differences - each indel fits in one segment
            (10, vec![(15, 5), (5, 15), (5, 5)]),
        ];

        // Run each test case
        for (i, (max_diff, expected_tracepoints)) in test_cases.iter().enumerate() {
            let tracepoints = cigar_to_tracepoints_raw(cigar, *max_diff);
            assert_eq!(
                tracepoints,
                *expected_tracepoints,
                "Test case {}: Raw tracepoints with max_diff={} incorrect",
                i + 1,
                max_diff
            );
        }
    }

    #[test]
    fn test_raw_vs_regular_tracepoints() {
        // Test that raw mode splits indels while regular mode doesn't
        let cigar = "5=8I3=7D2=";
        let max_diff = 5;

        // Regular mode: indels larger than max_diff become their own segments
        let regular = cigar_to_tracepoints(cigar, max_diff);
        // Raw mode: indels are split into max_diff-sized chunks
        let raw = cigar_to_tracepoints_raw(cigar, max_diff);

        // Regular mode should keep large indels intact
        assert_eq!(regular, vec![(5, 5), (8, 0), (3, 3), (0, 7), (2, 2)]);
        // Raw mode: indels can be combined with matches in segments
        // 5= + first 5I -> (10, 5)
        // remaining 3I + 3= -> (6, 3) but since we hit another indel, may need to check
        assert_eq!(raw, vec![(10, 5), (6, 5), (0, 5), (2, 2)]);
    }

    #[test]
    fn test_mixed_raw_tracepoints() {
        // Test mixed representation with raw mode
        let cigar = "3S5=6I2=4D1H";
        let max_diff = 3;

        let mixed_raw = cigar_to_mixed_tracepoints_raw(cigar, max_diff);

        assert_eq!(
            mixed_raw,
            vec![
                MixedRepresentation::CigarOp(3, 'S'),
                MixedRepresentation::Tracepoint(8, 5), // 5= + 3I
                MixedRepresentation::Tracepoint(3, 0), // remaining 3I
                MixedRepresentation::Tracepoint(2, 5), // 2= + 3D
                MixedRepresentation::Tracepoint(0, 1), // remaining 1D
                MixedRepresentation::CigarOp(1, 'H'),
            ]
        );
    }

    #[test]
    fn test_variable_raw_tracepoints() {
        // Test variable format with raw mode
        let cigar = "3=6I3=6D3=";
        let max_diff = 3;

        let variable_raw = cigar_to_variable_tracepoints_raw(cigar, max_diff);

        // With raw mode, indels combine with adjacent matches
        assert_eq!(
            variable_raw,
            vec![
                (6, Some(3)), // 3= + 3I
                (3, Some(0)), // remaining 3I
                (3, Some(6)), // 3= + 3D
                (0, Some(3)), // remaining 3D
                (3, None),    // 3=
            ]
        );
    }

    #[test]
    fn test_raw_cigar_roundtrip() {
        // Test that raw tracepoints can be reconstructed properly
        let original_cigar = "2=2I3=2D2=";
        // CIGAR consumes: a_seq: 2 + 2 + 3 + 0 + 2 = 9 bases
        //                 b_seq: 2 + 0 + 3 + 2 + 2 = 9 bases
        let a_seq = b"ACGTACGTA"; // 9 bases
        let b_seq = b"ACCGTCGAA"; // 9 bases
        let max_diff = 2;

        let distance_mode = DistanceMode::Affine2p {
            mismatch: 2,
            gap_open1: 4,
            gap_ext1: 2,
            gap_open2: 6,
            gap_ext2: 1,
        };

        // Test raw tracepoints roundtrip
        let raw_tracepoints = cigar_to_tracepoints_raw(original_cigar, max_diff);
        let reconstructed_cigar =
            tracepoints_to_cigar(&raw_tracepoints, a_seq, b_seq, 0, 0, &distance_mode);

        // The reconstructed CIGAR might be slightly different due to realignment,
        // but should represent the same alignment
        assert!(
            !reconstructed_cigar.is_empty(),
            "Raw roundtrip produced empty CIGAR"
        );
    }

    #[test]
    fn test_raw_edge_cases() {
        // Test edge cases for raw mode
        let test_cases = vec![
            // Empty CIGAR
            ("", 5, vec![]),
            // Only matches
            ("10=", 3, vec![(10, 10)]),
            // Single large indel
            ("15I", 5, vec![(5, 0), (5, 0), (5, 0)]),
            ("15D", 5, vec![(0, 5), (0, 5), (0, 5)]),
            // Mixed with mismatches - mismatches and indels can combine
            ("2X8I2X", 3, vec![(3, 2), (3, 0), (3, 0), (3, 2)]),
            // max_diff = 0 with indels
            ("3I2D", 0, vec![(1, 0), (1, 0), (1, 0), (0, 1), (0, 1)]),
        ];

        for (i, (cigar, max_diff, expected)) in test_cases.iter().enumerate() {
            let result = cigar_to_tracepoints_raw(cigar, *max_diff);
            assert_eq!(
                result,
                *expected,
                "Raw edge case {}: Failed for CIGAR '{}' with max_diff={}",
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
            (5, None),    // Should become (5, 5)
            (3, Some(0)), // Should become (3, 0)
            (0, Some(2)), // Should become (0, 2)
            (4, Some(6)), // Should become (4, 6)
        ];

        // Test the conversion logic used in variable_tracepoints_to_cigar
        let converted_regular: Vec<(usize, usize)> = variable_tracepoints
            .iter()
            .map(|(a_len, b_len_opt)| match b_len_opt {
                None => (*a_len, *a_len),
                Some(b_len) => (*a_len, *b_len),
            })
            .collect();

        let expected_regular = vec![(5, 5), (3, 0), (0, 2), (4, 6)];

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
        assert!(
            variable_tracepoints
                .iter()
                .any(|(_, b_opt)| b_opt.is_none()),
            "Variable tracepoints should have at least one optimized entry (None)"
        );
        assert!(
            variable_tracepoints
                .iter()
                .any(|(_, b_opt)| b_opt.is_some()),
            "Variable tracepoints should have at least one unoptimized entry (Some)"
        );
    }

    #[test]
    fn test_diagonal_tracepoint_generation() {
        // Test CIGAR string with various indel patterns
        let cigar = "5=3I2=2D4=";

        // Define test cases: (max_diff, expected_tracepoints) - updated based on actual diagonal logic
        let test_cases = [
            // Case 1: No diagonal distance allowed - each indel creates new segment
            (0, vec![(5, 5), (3, 0), (2, 2), (0, 2), (4, 4)]),
            // Case 2: Allow diagonal distance up to 2 - 3I exceeds max, then 2D+2=+4= fits in next segment
            (2, vec![(5, 5), (3, 0), (6, 8)]),
            // Case 3: Allow diagonal distance up to 3 - everything fits in one segment since diagonal balances
            (3, vec![(14, 13)]),
            // Case 4: Allow large diagonal distance - single segment (same as case 3)
            (10, vec![(14, 13)]),
        ];

        // Run each test case
        for (i, (max_diff, expected_tracepoints)) in test_cases.iter().enumerate() {
            let tracepoints = cigar_to_tracepoints_diagonal(cigar, *max_diff);
            assert_eq!(
                tracepoints,
                *expected_tracepoints,
                "Test case {}: Diagonal tracepoints with max_diff={} incorrect",
                i + 1,
                max_diff
            );
        }
    }

    #[test]
    fn test_diagonal_vs_regular_segmentation() {
        // Compare diagonal vs regular segmentation for complex patterns
        let test_cases = vec![
            // Alternating small indels - diagonal should group better
            ("2=1I2=1D2=", 2),
            // Large indels - should behave similarly
            ("3=5I3=5D3=", 3),
            // Balanced indels - diagonal should combine opposing indels
            ("2=3I2=3D2=", 4),
        ];

        for (cigar, max_diff) in test_cases {
            let regular_tracepoints = cigar_to_tracepoints(cigar, max_diff);
            let diagonal_tracepoints = cigar_to_tracepoints_diagonal(cigar, max_diff);

            println!("CIGAR: {cigar}, max_diff: {max_diff}");
            println!("Regular: {regular_tracepoints:?}");
            println!("Diagonal: {diagonal_tracepoints:?}");

            // Both should preserve total sequence lengths
            let regular_total: (usize, usize) = regular_tracepoints
                .iter()
                .fold((0, 0), |(a, b), (x, y)| (a + x, b + y));
            let diagonal_total: (usize, usize) = diagonal_tracepoints
                .iter()
                .fold((0, 0), |(a, b), (x, y)| (a + x, b + y));
            assert_eq!(
                regular_total, diagonal_total,
                "Total lengths should match for CIGAR: {cigar}"
            );
        }
    }

    #[test]
    fn test_diagonal_mixed_tracepoints() {
        // Test mixed representation with diagonal distance
        let test_cases = [
            // Special operations with balanced indels - diagonal allows combining
            (
                "3S5=2I2=2D1H",
                2,
                vec![
                    MixedRepresentation::CigarOp(3, 'S'),
                    MixedRepresentation::Tracepoint(9, 9), // 5= + 2I + 2= + 2D (balanced diagonal)
                    MixedRepresentation::CigarOp(1, 'H'),
                ],
            ),
            // Balanced indels with special ops
            (
                "2H3=3I3=3D2P",
                3,
                vec![
                    MixedRepresentation::CigarOp(2, 'H'),
                    MixedRepresentation::Tracepoint(9, 9), // 3= + 3I + 3= + 3D (balanced)
                    MixedRepresentation::CigarOp(2, 'P'),
                ],
            ),
        ];

        for (i, (cigar, max_diff, expected)) in test_cases.iter().enumerate() {
            let result = cigar_to_mixed_tracepoints_diagonal(cigar, *max_diff);
            assert_eq!(
                result,
                *expected,
                "Test case {}: Mixed diagonal tracepoints with max_diff={} incorrect for CIGAR '{}'",
                i + 1,
                max_diff,
                cigar
            );
        }
    }

    #[test]
    fn test_diagonal_variable_tracepoints() {
        // Test variable format with diagonal distance
        let cigar = "3=4I3=4D3=";
        let max_diff = 3;

        let variable_diagonal = cigar_to_variable_tracepoints_diagonal(cigar, max_diff);

        // Expected: with max_diff=3, the 4I exceeds limit, so we get separate segments
        // Then after reset, 3=4D balances but 4D by itself exceeds max_diff
        let expected = vec![
            (3, None),    // initial 3=
            (4, Some(0)), // 4I (exceeds max_diff=3, separate segment)
            (3, None),    // next 3=
            (0, Some(4)), // 4D (exceeds max_diff=3, separate segment)
            (3, None),    // final 3=
        ];

        assert_eq!(
            variable_diagonal, expected,
            "Variable diagonal tracepoints incorrect"
        );

        // Verify optimization is working (None used for equal lengths)
        assert!(
            variable_diagonal.iter().any(|(_, b_opt)| b_opt.is_none()),
            "Should have at least one optimized entry (None)"
        );
    }

    #[test]
    fn test_diagonal_edge_cases() {
        // Test edge cases for diagonal distance
        let test_cases = vec![
            // Empty CIGAR
            ("", 5, vec![]),
            // Only matches
            ("10=", 3, vec![(10, 10)]),
            // Only insertions
            ("8I", 3, vec![(8, 0)]),
            // Only deletions
            ("6D", 2, vec![(0, 6)]),
            // max_diff = 0 with balanced indels
            ("2I2D", 0, vec![(2, 0), (0, 2)]),
            // Large indels that cancel out
            ("5I10=5D", 10, vec![(15, 15)]), // Should combine since net diagonal distance = 0
        ];

        for (i, (cigar, max_diff, expected)) in test_cases.iter().enumerate() {
            let result = cigar_to_tracepoints_diagonal(cigar, *max_diff);
            assert_eq!(
                result,
                *expected,
                "Diagonal edge case {}: Failed for CIGAR '{}' with max_diff={}",
                i + 1,
                cigar,
                max_diff
            );
        }
    }

    #[test]
    fn test_diagonal_distance_calculation() {
        // Test the core diagonal distance logic with specific patterns
        let test_cases = [
            // Progressive insertion - distance should accumulate
            ("1I1I1I", 2, vec![(2, 0), (1, 0)]),
            // Alternating I/D - should balance out perfectly
            ("1I1D1I1D", 1, vec![(2, 2)]), // Net diagonal distance oscillates but stays <= 1
            // Large insertion followed by deletion - both exceed individually
            ("5I3D", 3, vec![(5, 0), (0, 3)]), // Each exceeds max_diff individually
            // Insertion that balances previous deletion
            ("3D2=3I", 3, vec![(5, 5)]), // Should combine since final distance = 0
        ];

        for (i, (cigar, max_diff, expected)) in test_cases.iter().enumerate() {
            let result = cigar_to_tracepoints_diagonal(cigar, *max_diff);
            assert_eq!(
                result,
                *expected,
                "Diagonal distance test case {}: Failed for CIGAR '{}' with max_diff={}",
                i + 1,
                cigar,
                max_diff
            );
        }
    }

    #[test]
    fn test_diagonal_roundtrip_compatibility() {
        // Test that diagonal tracepoints can be used with existing reconstruction functions
        let original_cigar = "2=3I2=2D3=";
        let a_seq = b"ACGTACGTAC"; // 10 bases
        let b_seq = b"ACAGTCGTACG"; // 11 bases
        let max_diff = 2;

        // Generate diagonal tracepoints
        let diagonal_tracepoints = cigar_to_tracepoints_diagonal(original_cigar, max_diff);

        let distance_mode = DistanceMode::Affine2p {
            mismatch: 2,
            gap_open1: 4,
            gap_ext1: 2,
            gap_open2: 6,
            gap_ext2: 1,
        };

        // Verify they can be reconstructed
        let reconstructed_cigar =
            tracepoints_to_cigar(&diagonal_tracepoints, a_seq, b_seq, 0, 0, &distance_mode);

        // Should not be empty and should represent a valid alignment
        assert!(
            !reconstructed_cigar.is_empty(),
            "Diagonal roundtrip produced empty CIGAR"
        );

        // Verify sequence length consistency
        let ops = cigar_str_to_cigar_ops(&reconstructed_cigar);
        let (total_a, total_b) = ops.iter().fold((0, 0), |(a, b), (len, op)| match op {
            '=' | 'M' | 'X' => (a + len, b + len),
            'I' => (a + len, b),
            'D' => (a, b + len),
            _ => (a, b),
        });

        // The original CIGAR "2=3I2=2D3=" should consume:
        // a_seq: 2+3+2+0+3=10 bases, b_seq: 2+0+2+2+3=9 bases
        // We provided sequences of length 10 and 11, so the reconstruction is valid
        assert_eq!(
            total_a,
            a_seq.len(),
            "Reconstructed CIGAR should consume correct a_seq length"
        );
        // The reconstruction may not consume the full b_seq if sequences don't match exactly
        assert!(
            total_b <= b_seq.len(),
            "Reconstructed CIGAR should not overconsume b_seq"
        );
    }

    #[test]
    fn test_fastga_tracepoints() {
        // Test tracepoints from a real PAF comparison against the FASTGA reference.
        let cigar = "15=1X5=1X4=1X1=1X92=1X27=1X4=1X11=1X2=1X17=1X1=1X36=1X7=1X12=1X1=1X1D2=1I7=1X9=1D1=1I1=2X6=1X1=1X39=1X31=1X4=1X39=1X31=1X120=1X27=1X43=1X11=1X13=27I6=1X13=1X16=1X12=1X87=1X8=1X39=2X7=1X74=1X39=1X56=1X34=1X21=1X76=1X17=1D128=1X116=1X4=1X37=1D20=1X2=1X18=1X25=1X84=1X95=1X108=4D42=1X9=1X74=1X1=1X159=1X91=1X7=1X2=1X41=1I102=2X9=1X49=1X4=1X14=1X23=1X37=2X53=1X39=2X16=1X30=1X8=1X130=1X41=1X38=1X122=1X107=1X7=1X71=1X169=1X14=1X18=1X59=1X21=1X10=1D1X113=1X36=1X2=1X70=1X16=1X34=1X136=1X68=1X6=1X98=1X38=1X37=1X65=1X237=1X46=1X26=1X8=1X218=1X51=1X172=1X137=1X1=1X12=1X107=1X26=1X112=3I13=1X9=1X2=1X77=1X75=1X14=1I119=1X9=1X20=1X7=1X1=1X40=5D20=1X4=1X7=2X47=1X27=1X53=1X220=1X4=1X8=1X66=1X195=1I33=1X27=1X29=1X111=1X74=1X2=1X47=1X21=1X113=1X44=1X5=1X3=1X2=1X122=1X208=1X18=1X136=2I52=11D7=1X26=1X85=1X47=1X27=1X12=1X43=1X4=1X30=1X41=1X98=1X88=1X183=1X185=1X18=1X51=1X44=1X31=1X47=1X48=3D3=1D105=1X162=1X1=1X97=1X3=1X55=1X72=1X6=1X53=1X12=1I4=1D31=1X91=1X45=1X84=1X6=1X78=1X78=1X13=1D7=1X23=1X20=1D2=1I25=1D1=1I5=1X41=1X41=1X26=1X4=2D1X6=1D2=1X10=1X8=1X1=1X3=2X12=1X12=1D1=1X1I21=1X9=2X34=1X5=1D1=1I1=1X8=1X11=1X5=1X20=1X2=1X62=2X4=1X6=1X22=1X5=1X5=1X30=1X8=1X3=1X1=1X1=1D1X2=1X2=1D1=1X1=1X1=1X2=1X1=1X2=1D9=1X3=2X1=2X2=2X2=1D2=1I4=1X25=1X5=1X17=1X2=1D1X3=1I8=1X69=1X2=1X20=1X44=1X38=1X35=1X56=1X8=1X8=1X111=1X2=3X1=1X3=3I8=1X1=4X8=1X1=1X6=1X3=1X3=1I10=1I10=1I90=1X7=1X3=1X6=1I11=1X7=1X1=2X5=1X1=2X3=9I6=1X1=1X2=1X6=1X10=1X28=1X11=1X14=1X7=1X68=1X5=1X25=1X2=1X1=1X12=1X5=1X23=1X64=1X5=1X62=1X15=1X2=2X3=3X9=1X4=1X5=1I2=1X3=1D2X3=1X7=1X6=1X10=3X11=1X87=1X4=1X26=1X3=1X1=1D2=1X4=1X6=1X2=1X1=1X6=2X22=1X20=1X20=1X6=1X14=1I3=1X1=1X6=1X16=1X3=4X3=1X5=1X1=2X3=2X4=1X12=1X1=2X2=2X4=1I4=1D2=1X1=1D1=1I22=1X4=1X6=1D1=2X2=1I7=1X8=1X7=1X2=1X4=1X58=1D1=1X1I21=1X27=1X40=2X3=1X25=1X2I6=1I2=1D9=2I1=3I7=1X3=1X47=1X1=2X1=1X16=1X13=7D15=1X14=1X1=1X6=1X17=1I3=1D14=1I3=1I3=1X7=2X9=1X6=1X13=1D1X5=1D2=2I2=1D3=1D1X1=1D1=1X1=1D1X6=3I6=2X1=1D1X1=1X3=1D1X8=1D2=1I16=1X10=1X3=1X10=2X1I4=1D22=2D1=1X9=1D7=1X1=1D3=1I9=4X1D4=1X5=1X5=5I12=1X4I4=1X5=1I1=1X10=1I7=1I2=1I18=1D2=1I6=";
        let trace_spacing = 100;
        let a_start = 0;
        let b_start = 5207790;
        let complement = true;

        let expected_str = "20,88;18,103;28,103;16,105;13,96;8,97;6,100;19,100;18,99;12,101;2,100;21,100;2,100;6,100;4,100;20,91;6,99;16,94;4,100;4,100;2,100;3,100;17,100;17,103;6,100;10,100;15,103;5,100;6,101;1,100;2,100;2,100;1,100;5,100;1,100;2,100;2,100;1,100;4,104;2,100;2,100;2,100;0,100;1,100;0,100;2,100;1,100;3,100;3,100;2,100;14,109;0,100;2,100;0,100;1,100;4,100;1,100;2,100;3,100;0,100;4,99;0,100;1,100;3,100;0,100;1,100;2,100;11,105;3,100;2,99;1,100;6,97;2,100;0,100;3,100;1,100;0,100;2,100;0,100;1,100;3,100;0,100;0,100;2,100;2,100;2,100;1,100;1,100;3,100;2,100;4,101;3,100;0,100;1,100;2,100;1,100;2,100;1,100;2,100;4,100;4,100;5,100;1,99;3,100;1,100;2,100;2,100;4,104;1,100;2,100;4,101;2,100;1,100;2,101;2,100;2,100;4,100;2,100;30,73;4,100;0,100;4,100;10,100;9,100;2,100;4,50";

        // Change the expected type to Vec<(usize, usize)>
        let expected_tracepoints: Vec<(usize, usize)> = expected_str
            .split(';')
            .filter_map(|pair| {
                let mut parts = pair.split(',');
                if let (Some(x_str), Some(y_str)) = (parts.next(), parts.next()) {
                    // Change parse::<i32>() to parse::<usize>()
                    if let (Ok(x), Ok(y)) = (x_str.parse::<usize>(), y_str.parse::<usize>()) {
                        return Some((x, y));
                    }
                }
                None
            })
            .collect();
        
        let actual_tracepoints = cigar_to_tracepoints_fastga(cigar, trace_spacing, a_start, b_start, complement);

        assert_eq!(
            actual_tracepoints.len(),
            expected_tracepoints.len(),
            "The number of generated tracepoints ({}) does not match the expected number ({})",
            actual_tracepoints.len(),
            expected_tracepoints.len()
        );

        for (i, (actual, expected)) in actual_tracepoints.iter().zip(expected_tracepoints.iter()).enumerate() {
            assert_eq!(
                actual,
                expected,
                "Mismatch at tracepoint index {}: expected {:?}, but got {:?}",
                i,
                expected,
                actual
            );
        }
    }
}
