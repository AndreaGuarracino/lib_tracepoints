# lib_tracepoints
A Rust library for compressed alignment representation using tracepoints.

## Overview

`lib_tracepoints` provides utilities for converting between CIGAR strings and variable-delta tracepoints for efficient alignment representation in bioinformatics applications. The library enables:

- Converting CIGAR strings to variable-delta tracepoints
- Reconstructing CIGAR strings from tracepoints
- Managing alignment segments with configurable difference thresholds

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

```rust
use lib_tracepoints::{cigar_to_tracepoints_variable, tracepoints_to_cigar_variable};

fn main() {
    // Convert CIGAR to tracepoints with max difference of 5
    let cigar = "10=2D5=2I3=";
    let tracepoints = cigar_to_tracepoints_variable(&cigar, 5);
    
    // Reconstruct CIGAR from tracepoints
    let a_seq = "ACGTACGTACACGTACGTAC";
    let b_seq = "ACGTACGTACACGTACGTAC";
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
}
```