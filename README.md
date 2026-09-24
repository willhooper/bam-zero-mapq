# bam-zero-unmapped-mapq

A fast streaming Rust CLI that sets `MAPQ` to `0` whenever the BAM unmapped
flag (`0x4`) is set. Every other field and the input record order are preserved.

It uses HTSlib's native shared thread pool to overlap BGZF decompression and
compression. It builds a BAI index during the same write pass, using HTSlib's
`sam_idx_init` / `sam_write1` / `sam_idx_save` APIs. Index entries are accumulated
as compressed output blocks are written, then saved to `<output>.bai` at completion.
There is no second pass over the BAM. The record loop reuses one allocation and
does not load the BAM into memory.

## Build

You need a Rust toolchain and a C compiler. Then run:

```sh
cargo build --release
```

The dependency is intentionally pinned to a rust-htslib release that uses its
supplied Linux/macOS bindings, so building does not require libclang/bindgen.
Its matching `hts-sys` version is pinned too because newer generated bindings
are not source-compatible with that rust-htslib release.

## Run

```sh
./target/release/bam-zero-unmapped-mapq \
  --threads 16 input.bam output.bam
```

This creates both `output.bam` and `output.bam.bai`. Input must already be
coordinate-sorted, with reads without coordinates at the end. BAI supports
reference coordinates below 2^29 (512 Mbp). Unsorted or out-of-range records,
index failures, and output write failures cause a nonzero exit status.

To sort upstream in a pipeline:

```sh
samtools sort -@ 16 -O BAM input.bam \
  | ./target/release/bam-zero-unmapped-mapq -t 16 - output.bam
```

Input may be `-` for stdin. Output must be a file path; stdout output is no
longer supported so every successful run produces a BAM and its BAI index.

Output compression defaults to level 1 for throughput. Use `--compression 0`
for the highest throughput and larger output, or a value up to 9 for a smaller
file at the cost of speed.

Run `bam-zero-unmapped-mapq --help` for all options.
