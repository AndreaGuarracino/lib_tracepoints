# tracepoints

A library for sequence alignment compression and reconstruction using tracepoints.

## Overview

`tracepoints` provides utilities for converting between CIGAR strings and tracepoint representations for efficient alignment storage. The library enables:

- Converting CIGAR strings to tracepoints
- Reconstructing CIGAR strings from tracepoints

### What are tracepoints?

Rather than storing every alignment operation in a full CIGAR string, tracepoints record a sparse set of coordinate pairs along the alignment path. Each pair of consecutive tracepoints defines a short subalignment interval whose CIGAR can be reconstructed on-demand by re-aligning the corresponding sequence segments with WFA.

This library implements **adaptive tracepoints**: instead of segmenting at fixed intervals, it segments based on local alignment complexity, creating larger segments in conserved regions and smaller ones in divergent regions. Reconstruction from adaptive tracepoints guarantees identical or improved alignment scores, never worse.

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

- **Tracepoint types**:
  - **Standard**: `(a_len, b_len)` pairs for each segment
  - **FastGA**: Fixed-spacing tracepoints compatible with the FastGA aligner
- **Complexity metrics**: `EditDistance` (count of mismatches + indels) and `DiagonalDistance` (max diagonal shift within a segment)
- **CIGAR reconstruction**: Conversion from tracepoints back to CIGAR strings using WFA alignment
- **Distance modes**: Support for edit distance, gap-affine, and dual gap-affine penalties

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
