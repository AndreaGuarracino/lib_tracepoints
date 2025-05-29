# lib_tracepoints

A Rust library for sequence alignment compression and reconstruction using tracepoints.

## Overview

`lib_tracepoints` provides utilities for converting between CIGAR strings and various tracepoint representations for efficient alignment storage in bioinformatics applications. The library enables:

- Converting CIGAR strings to tracepoints
- Reconstructing CIGAR strings from tracepoints
- Supporting multiple tracepoint representations with varying space-time tradeoffs

## Installation

Add this to the `Cargo.toml` file of the project where you want to use `lib_tracepoints`:

```toml
[dependencies]
lib_wfa2 = { git = "https://github.com/AndreaGuarracino/lib_wfa2"}
lib_tracepoints = { git = "https://github.com/AndreaGuarracino/lib_tracepoints"}
```

This library depends on `lib_wfa2`, which requires the `WFA2-lib` (commit `3c1734e9bb319c7782ae6845e627612ff157d1cc`) to be built first:

```shell
git clone https://github.com/smarco/WFA2-lib
cd WFA2-lib
git checkout 3c1734e9bb319c7782ae6845e627612ff157d1cc
make clean all
```

Then build the project where you want to use `lib_tracepoints` with:

```shell
# Point to your pre-built WFA2-lib directory
export WFA2LIB_PATH="/path/to/WFA2-lib"
# Build your project
cargo build --release
```

## Features

- **Basic Tracepoints**: Simple `(a_len, b_len)` pairs for each segment
- **Double-Band Tracepoints**: Enhanced `(a_len, b_len, (min_diagonal, max_diagonal))` triples for faster reconstruction
- **Single-Band Tracepoints**: Memory-efficient `(a_len, b_len, max_abs_diagonal)` triples
- **Variable-Band Tracepoints**: Adaptive representation that optimizes for different diagonal patterns
- **Mixed Representation**: Preserves special CIGAR operations `(H, N, P, S)` that aren't suitable for reconstruction

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

### Single-band tracepoints

This approach stores only the maximum absolute diagonal boundary for each tracepoint. It's more memory-efficient but may require more computation during alignment reconstruction.

```rust
use lib_tracepoints::{cigar_to_single_band_tracepoints, single_band_tracepoints_to_cigar};

fn main() {
    let cigar = "10=2D5=2I3=";
    
    // Store only the maximum absolute diagonal value (more compact)
    let single_band_tracepoints = cigar_to_single_band_tracepoints(cigar, 5);
    
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

### Double-band tracepoints

This approach stores both the minimum and maximum diagonal boundaries for each tracepoint. It provides more detailed boundary information for faster CIGAR string reconstruction but results in a slightly larger alignment representation.

```rust
use lib_tracepoints::{cigar_to_double_band_tracepoints, double_band_tracepoints_to_cigar};

fn main() {
    let cigar = "10=2D5=2I3=";
    
    // Store both min and max diagonal boundaries
    let tracepoints = cigar_to_double_band_tracepoints(cigar, 5);

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

### Variable-band tracepoints

This approach adaptively selects the most efficient representation based on diagonal patterns, optimizing storage while maintaining alignment reconstruction speed.

```rust
use lib_tracepoints::{cigar_to_variable_band_tracepoints, variable_band_tracepoints_to_cigar};

fn main() {
    let cigar = "10=2D5=2I3=";
    
    // Use an optimized representation based on diagonal patterns
    let variable_band_tracepoints = cigar_to_variable_band_tracepoints(cigar, 5);
    
    let a_seq = "ACGTACGTACACGTACGTAC";
    let b_seq = "ACGTACGTACACGTACGTAC";
    let reconstructed_cigar = variable_band_tracepoints_to_cigar(
        &variable_band_tracepoints,
        a_seq.as_bytes(),
        b_seq.as_bytes(),
        0, 0, 2, 4, 2, 6, 1
    );
}
```

### Mixed representation

This approach preserves special CIGAR operations that aren't suitable for WFA alignment.

```rust
use lib_tracepoints::{cigar_to_mixed_double_band_tracepoints, mixed_double_band_tracepoints_to_cigar, MixedRepresentation};

fn main() {
    let cigar = "5=2H3=";
    let mixed_tracepoints = cigar_to_mixed_double_band_tracepoints(cigar, 2);
    
    let a_seq = b"ACGTACGTAC";
    let b_seq = b"ACGTACGTAC";
    let reconstructed_cigar = mixed_double_band_tracepoints_to_cigar(
        &mixed_tracepoints,
        a_seq,
        b_seq,
        0, 0, 2, 4, 2, 6, 1
    );
}
```

## How It Works

### CIGAR to Tracepoints Conversion

The library segments a CIGAR string into tracepoints where each segment contains at most `max_diff` differences (mismatches or indels):

- Match operations ('=' and 'M') don't count as differences
- Mismatch operations ('X') can be split across segments if needed
- Indels ('I', 'D') are kept intact within a single segment when possible
- Long indels exceeding `max_diff` become their own segments

### Diagonal Tracking

For banded implementations:
- `min_diagonal` tracks the lowest diagonal reached in a segment
- `max_diagonal` tracks the highest diagonal reached in a segment
- Diagonal position changes when indels are encountered (increases for insertions, decreases for deletions)

### Tracepoints to CIGAR Reconstruction

For each tracepoint pair, the library performs the alignment of the corresponding sequence segments using WFA alignment:
- Pure insertions (a_len > 0, b_len = 0) are directly converted to 'I' operations
- Pure deletions (a_len = 0, b_len > 0) are directly converted to 'D' operations
- Mixed segments are realigned using the WFA algorithm with the diagonal boundary information
