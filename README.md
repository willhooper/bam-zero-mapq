# bam-zero-unmapped-mapq

A fast streaming Rust CLI that sets `MAPQ` to `0` whenever the BAM unmapped
flag (`0x4`) is set. Every other field and the input record order are preserved.

It uses HTSlib's native shared thread pool to overlap BGZF decompression and
compression. The record loop reuses one allocation and does not load the BAM
into memory.

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

For a pipeline:

```sh
samtools view -u input.bam \
  | ./target/release/bam-zero-unmapped-mapq -t 16 - - \
  | samtools sort -@ 16 -o output.bam -
```

Output compression defaults to level 1 for throughput. Use `--compression 0`
for the highest throughput and larger output, or a value up to 9 for a smaller
file at the cost of speed.

Run `bam-zero-unmapped-mapq --help` for all options.
