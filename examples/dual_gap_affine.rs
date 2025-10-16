use lib_tracepoints::{cigar_to_tracepoints, tracepoints_to_cigar, DistanceMode};

fn main() {
    let a_seq = b"ACGTACGTACGT"; // 12 bases
    let b_seq = b"ATGTACCTAC";   // 10 bases
    let cigar = "1=1X4=1X3=2I";

    // Convert to tracepoints with max 2 differences per segment
    let tracepoints = cigar_to_tracepoints(cigar, 2);

    // Use affine gap penalties
    let distance_mode = DistanceMode::Affine2p {
        mismatch: 2,
        gap_open1: 4,
        gap_ext1: 2,
        gap_open2: 6,
        gap_ext2: 1,
    };

    // Reconstruct CIGAR
    let reconstructed = tracepoints_to_cigar(&tracepoints, a_seq, b_seq, 0, 0, &distance_mode);

    println!("Tracepoints: {:?}", tracepoints);
    println!("     Original CIGAR: {}", cigar);
    println!("Reconstructed CIGAR: {}", reconstructed);
}
