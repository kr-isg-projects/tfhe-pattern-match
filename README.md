# Efficient Homomorphic String Search via TFHE

This repository provides the source code used in the paper *"Efficient Homomorphic String Search via TFHE"*.

## Requirements
- SageMath (version 9.5 or later)
- Rust (rustc version 1.93 or later, for executing the TFHE circuits)

## Dependencies
Our Rust code depends on:  
- [TFHE-rs v0.5.3](https://github.com/zama-ai/tfhe-rs/tree/tfhe-rs-0.5.3)  
- [refined-tfhe-lhe](https://github.com/KAIST-CryptLab/refined-tfhe-lhe) (commit 9b0426e) with some modifications  

Please clone or download these repositories and place them in `tfhe-rs` and `refined-tfhe-lhe` directories respectively. Then, replace or patch the following files as indicated.

### Modifications
- In **TFHE-rs v0.5.3**:  
  `tfhe-rs/tfhe/src/core_crypto/fft_impl/mod.rs` — additional functions added  

- In **refined-tfhe-lhe**:  
  `refined-tfhe-lhe/error_analysis/fhe_pattern.sage` — modifications based on `integer_input_lhe.sage`  

---

## Reproducing Experimental Results
To reproduce results shown in Table 2 and Table 3 of the paper, execute:

```bash
/usr/bin/time -v cargo run --release -- \
    --seed 0 --log-n 12 --n0 1080 --m 100 --m0 100 --solution 0
```
Where:

--log-n : text length in log₂

--n0 : text length without padding

--m : pattern length

--m0 : pattern length without padding


## TFHE Parameter Search

To search for TFHE parameters suitable for our circuit, execute:

```bash
sage fhe_pattern.sage
```

## License

This software is distributed under the BSD-3-Clause-Clear license.

## Patents and Commercial Use

This software may be covered by one or more patents.
For commercial use of this software or any underlying patented technology,
a separate patent license may be required.

For inquiries regarding commercial use, please contact:
sh-narisada.at.kddi.com

## Funding
This work was partially supported by JST K Program, grant number JPMJKP24U2, Japan.