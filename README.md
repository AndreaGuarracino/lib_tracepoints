# lib_tracepoints

A Rust library for compressed alignment representation using tracepoints.

## Overview

`lib_tracepoints` provides utilities for converting between CIGAR strings and tracepoints for efficient alignment representation in bioinformatics applications. The library enables:

- Converting CIGAR strings to tracepoints
- Reconstructing CIGAR strings from tracepoints
- Managing alignment segments with configurable difference thresholds
- Tuning alignment compression/decompression with diagonal boundary tracking (dual or symmetric)

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

### Banded tracepoints

The library supports two approaches for diagonal boundary tracking.

The full diagonal tracking stores both the minimum and maximum diagonal boundaries for each tracepoint. This approach is more verbose (compress less) but allows for faster alignment reconstruction.

```rust
use lib_tracepoints::{cigar_to_banded_tracepoints, banded_tracepoints_to_cigar};

fn main() {
    let cigar = "10=2D5=2I3=";
    
    // Store both min and max diagonal boundaries
    let tracepoints: Vec<(usize, usize, (isize, isize))> = 
        cigar_to_banded_tracepoints(cigar, 5);

    let a_seq = "ACGTACGTACACGTACGTAC";
    let b_seq = "ACGTACGTACACGTACGTAC";
    let reconstructed_cigar = banded_tracepoints_to_cigar(
        &tracepoints,
        a_seq.as_bytes(),
        b_seq.as_bytes(),
        0, 0, 2, 4, 2, 6, 1
    );
}
```

The symmetric band tracking stores only the maximum diagonal boundary for each tracepoint. This approach is more compact (compress more) but requires more computation to reconstruct the alignment compared to the full diagonal tracking.

```rust
use lib_tracepoints::{cigar_to_symmetric_banded_tracepoints, symmetric_banded_tracepoints_to_cigar};

fn main() {
    let cigar = "10=2D5=2I3=";
    
    // Store only the maximum absolute diagonal value (more compact)
    let symmetric_tracepoints: Vec<(usize, usize, usize)> = 
        cigar_to_symmetric_banded_tracepoints(cigar, 5);
    
    let a_seq = "ACGTACGTACACGTACGTAC";
    let b_seq = "ACGTACGTACACGTACGTAC";
    let reconstructed_cigar = symmetric_banded_tracepoints_to_cigar(
        &symmetric_tracepoints,
        a_seq.as_bytes(),
        b_seq.as_bytes(),
        0, 0, 2, 4, 2, 6, 1
    );
}
```
