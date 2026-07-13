# rsomics-fasta-translate

Translate a DNA/RNA FASTA to protein in one or more reading frames — a Rust
port of `seqkit translate`.

## Install

```
cargo install rsomics-fasta-translate
```

## Usage

```
# reading frame 1 (default)
rsomics-fasta-translate input.fasta

# all six frames
rsomics-fasta-translate input.fasta -f 6

# a chosen set of frames, to a file
rsomics-fasta-translate input.fasta -f 1,3,-2 -o proteins.faa
```

- `-f, --frames` — comma-separated frames (`1,2,3,-1,-2,-3`), or `6` for all
  six (default `1`).
- `-o, --output` — output path (`-` = stdout).

## Origin

Independent Rust reimplementation of `seqkit translate`, based on documented
behaviour and black-box comparison against `seqkit` v2.13, using the standard
genetic code with IUPAC-ambiguity-aware codon resolution. Frozen goldens live
in `tests/golden/`.

License: MIT OR Apache-2.0.
Upstream credit: [seqkit](https://github.com/shenwei356/seqkit) (MIT).
