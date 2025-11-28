# tracepoints

A library for sequence alignment compression and reconstruction using tracepoints.

## Overview

`tracepoints` provides utilities for converting between CIGAR strings and tracepoint representations for efficient alignment storage. The library enables:

- Converting CIGAR strings to tracepoints
- Reconstructing CIGAR strings from tracepoints

## Installation

Add this to your `Cargo.toml`:

```toml
[dependencies]
tracepoints = { git = "https://github.com/AndreaGuarracino/tracepoints" }
```

Then simply build your project:

```shell
cargo build --release
```

## Features

- **Tracepoints**: Simple `(a_len, b_len)` pairs for each segment
- **CIGAR reconstruction**: Efficient conversion from tracepoints back to CIGAR strings using WFA alignment
- **Distance modes**: Support for both affine gap penalties and edit distance

## Usage

See the examples:
- [`dual_gap_affine.rs`](examples/dual_gap_affine.rs) - Using dual gap-affine distance
- [`edit_distance.rs`](examples/edit_distance.rs) - Using edit distance

Run with: `cargo run --example dual_gap_affine` or `cargo run --example edit_distance`.

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
