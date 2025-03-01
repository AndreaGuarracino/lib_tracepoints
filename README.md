# lib_tracepoints

A Rust library for compressed alignment representation using tracepoints.

## Overview

`lib_tracepoints` provides utilities for converting between CIGAR strings and tracepoints for efficient alignment representation in bioinformatics applications. The library enables:

- Converting CIGAR strings to tracepoints
- Reconstructing CIGAR strings from tracepoints
- Managing tracepoints with configurable difference thresholds
- Tuning alignment compression/decompression with diagonal boundary tracking (single- or double-band)

## Installation

Add this to your `Cargo.toml`:

```toml
[dependencies]
lib_tracepoints = { git = "https://github.com/AndreaGuarracino/lib_tracepoints" }
```

## Dependencies

This library depends on `lib_wfa2`, which requires the `WFA2-lib` to be built first:

```shell
git clone https://github.com/smarco/WFA2-lib
cd WFA2-lib
make clean all
```

Then build your project with:

```shell
# Point to your pre-built WFA2-lib directory
export WFA2LIB_PATH="/path/to/WFA2-lib"
# Build your project
cargo build --release
```

## Usage

### Basic tracepoints

```rust
use lib_tracepoints::{cigar_to_tracepoints, tracepoints_to_cigar};

fn main() {
    // Convert CIGAR to tracepoints with max difference of 5
    let cigar = "10=2D5=2I3=";
    let tracepoints = cigar_to_tracepoints(&cigar, 5);
    
    // Reconstruct CIGAR from tracepoints
    let a_seq = "ACGTACGTACACGTACGTAC";
    let b_seq = "ACGTACGTACACGTACGTAC";
    let reconstructed_cigar = tracepoints_to_cigar(
        &tracepoints,
        a_seq.as_bytes(),
        b_seq.as_bytes(),
        0,  // a_start
        0,  // b_start
        2,  // mismatch penalty
        4,  // gap_open1
        2,  // gap_ext1
        6,  // gap_open2
        1,  // gap_ext2
    );
}
```

### Double-band tracepoints

This approach stores both the minimum and maximum diagonal boundaries for each tracepoint. It provides more detailed boundary information for faster CIGAR string reconstruction but results in a slightly larger alignment representation.

```rust
use lib_tracepoints::{cigar_to_double_band_tracepoints, double_band_tracepoints_to_cigar};

fn main() {
    let cigar = "10=2D5=2I3=";
    
    // Store both min and max diagonal boundaries
    let tracepoints: Vec<(usize, usize, (isize, isize))> = 
        cigar_to_double_band_tracepoints(cigar, 5);

    let a_seq = "ACGTACGTACACGTACGTAC";
    let b_seq = "ACGTACGTACACGTACGTAC";
    let reconstructed_cigar = double_band_tracepoints_to_cigar(
        &tracepoints,
        a_seq.as_bytes(),
        b_seq.as_bytes(),
        0, 0, 2, 4, 2, 6, 1
    );
}
```

### Single-band tracepoints

This approach stores only the maximum absolute diagonal boundary for each tracepoint. It's more memory-efficient but may require more computation during alignment reconstruction.

```rust
use lib_tracepoints::{cigar_to_single_band_tracepoints, single_band_tracepoints_to_cigar};

fn main() {
    let cigar = "10=2D5=2I3=";
    
    // Store only the maximum absolute diagonal value (more compact)
    let single_band_tracepoints: Vec<(usize, usize, usize)> = 
        cigar_to_single_band_tracepoints(cigar, 5);
    
    let a_seq = "ACGTACGTACACGTACGTAC";
    let b_seq = "ACGTACGTACACGTACGTAC";
    let reconstructed_cigar = single_band_tracepoints_to_cigar(
        &single_band_tracepoints,
        a_seq.as_bytes(),
        b_seq.as_bytes(),
        0, 0, 2, 4, 2, 6, 1
    );
}
```
