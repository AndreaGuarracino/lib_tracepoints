# lib_tracepoints

A library for sequence alignment compression and reconstruction using tracepoints.

## Overview

`lib_tracepoints` provides utilities for converting between CIGAR strings and tracepoint representations for efficient alignment storage. The library enables:

- Converting CIGAR strings to tracepoints
- Reconstructing CIGAR strings from tracepoints

## Installation

Add this to your `Cargo.toml`:

```toml
[dependencies]
lib_tracepoints = { git = "https://github.com/AndreaGuarracino/lib_tracepoints" }
```

Then simply build your project:

```shell
cargo build --release
```

## Features

- **Tracepoints**: Simple `(a_len, b_len)` pairs for each segment
- **CIGAR Reconstruction**: Efficient conversion from tracepoints back to CIGAR strings using WFA alignment

## Usage

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

## How It Works

### CIGAR to Tracepoints Conversion

The library segments a CIGAR string into tracepoints where each segment contains at most `max_diff` differences (mismatches or indels):

- Match operations ('=' and 'M') don't count as differences
- Mismatch operations ('X') can be split across segments if needed
- Indels ('I', 'D') are kept intact within a single segment when possible
- Long indels exceeding `max_diff` become their own segments


### Tracepoints to CIGAR Reconstruction

For each tracepoint pair, the library performs the alignment of the corresponding sequence segments using WFA alignment:
- Pure insertions (a_len > 0, b_len = 0) are directly converted to 'I' operations
- Pure deletions (a_len = 0, b_len > 0) are directly converted to 'D' operations
- Mixed segments are realigned using the WFA algorithm
