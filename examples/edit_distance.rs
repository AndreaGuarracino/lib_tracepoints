use lib_tracepoints::{cigar_to_tracepoints, tracepoints_to_cigar, DistanceMode};

fn main() {
    let cigar = "5=1I5=";
    let a_seq = b"ACGTACGTACG"; // 11 bases
    let b_seq = b"ACGTAGTACG"; // 10 bases

    // Convert to tracepoints with max 3 differences per segment
    let tracepoints = cigar_to_tracepoints(cigar, 3);

    // Use edit distance mode (unit costs)
    let edit_mode = DistanceMode::Edit {
        mismatch: 1,
        gap_opening: 1,
    };

    // Reconstruct CIGAR
    let reconstructed = tracepoints_to_cigar(&tracepoints, a_seq, b_seq, 0, 0, &edit_mode);

    println!("Tracepoints: {:?}", tracepoints);
    println!("     Original CIGAR: {}", cigar);
    println!("Reconstructed CIGAR: {}", reconstructed);
}
