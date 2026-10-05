# Efficient Homomorphic String Search via TFHE

This repository provides an implementation of the method proposed in the paper **"Efficient Homomorphic String Search via TFHE"**.

Paper: https://www.sciltp.com/journals/pc/articles/2608004913

## Features

- Search over encrypted text using TFHE
- `O(m log n)` secure character comparisons
- No leakage of the search result during server-side computation
- ASCII character support
- DNA sequence support (`A`, `C`, `G`, `T`)
- Arbitrary text length and pattern length

## Requirements

- Rust 1.93 or later
- SageMath 9.5 or later (for parameter analysis)

## Dependencies

- [TFHE-rs v1.8.1](https://github.com/zama-ai/tfhe-rs)

## Usage

### Search a pattern in ASCII text

ASCII is the default alphabet.

```bash
cargo run --release -- \
  --text 'Call me Ishmael. Some years ago - never mind how long precisely - having little or no money in my purse, and nothing particular to interest me on shore, I thought I would sail about a little and see the watery part of the world.' \
  --pattern 'having little or no money in my purse'
```

### Search a pattern in a DNA sequence

```bash
cargo run --release -- \
  --alphabet dna \
  --text 'TAGGACCTGATCGTACGATCG' \
  --pattern 'GACCTG'
```

### Search a pattern in a large DNA sequence

```bash
FHE_FIND_PROFILE=1 /usr/bin/time -v cargo run --release -- \
  --alphabet dna \
  --seed 0 \
  --log-n 20 \
  --n0 100000 \
  --m 100 \
  --m0 100
```

### Search a pattern from files

Example files are provided in `example/`.

```bash
cargo run --release -- \
  --text-file example/text.txt \
  --pattern-file example/pattern.txt
```

## TFHE Parameter Analysis

```bash
cd error-analysis
sage tfhe_pattern.sage
```

## Version History

See [Releases](https://github.com/kr-isg-projects/tfhe-pattern-match/releases).

### v0.1.0

- DNA sequence support (`A`, `C`, `G`, `T`)
- Used to reproduce the results reported in the paper
- 128-bit security with failure probability (FP) below `2^-60`
- Pattern length limited to `m <= N/2`
- Based on TFHE-rs v0.5.3 and `refined-tfhe-lhe`

### v0.2.0

- Added ASCII character support
- Added arbitrary pattern-length support
- 128-bit security with FP below `2^-128`
- Updated to TFHE-rs v1.8.1
- Removed the dependency on `refined-tfhe-lhe`
- Added direct and file-based text/pattern input

## License

This software is distributed under the BSD-3-Clause-Clear license.

## Patents and Commercial Use

This software may be covered by one or more patents.

For commercial use of this software or any underlying patented technology, a separate patent license may be required.

For inquiries regarding commercial use, please contact:

`sh-narisada.at.kddi.com`

## Funding

This work was partially supported by JST K Program, grant number `JPMJKP24U2`, Japan.