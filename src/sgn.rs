use rayon::prelude::*;
use tfhe::core_crypto::entities::{FourierLweBootstrapKey, LweCiphertext, LweSecretKey};
use tfhe::core_crypto::prelude::*;

use crate::pbs::*;

/// Return `right` if `left == 0`, otherwise return `left`.
///
/// Inputs and output are in {0, 1, 15}, where 15 represents -1 mod 16.
fn combine_first_nonzero(
    left: &LweCiphertext<Vec<u64>>,
    right: &LweCiphertext<Vec<u64>>,
    ciphertext_modulus: CiphertextModulus<u64>,
    lwe_secret_key_after_ks: &LweSecretKey<Vec<u64>>,
    ksk: &LweKeyswitchKey<Vec<u64>>,
    bsk: &FourierLweBootstrapKey<aligned_vec::ABox<[tfhe::core_crypto::fft_impl::fft64::c64]>>,
    polynomial_size: PolynomialSize,
    glwe_size: GlweSize,
    delta: u64,
) -> LweCiphertext<Vec<u64>> {
    // Encode the pair as:
    //
    // x = left + 3 * right
    //
    // Required mapping:
    //
    // left=0:
    //   right=0   -> x=0  -> 0
    //   right=1   -> x=3  -> 1
    //   right=15  -> x=13 -> 15
    //
    // left=1:
    //   right=0   -> x=1  -> 1
    //   right=1   -> x=4  -> 1
    //   right=15  -> x=14 -> 1
    //
    // left=15:
    //   right=0   -> x=15 -> 15
    //   right=1   -> x=2  -> 15
    //   right=15  -> x=12 -> 15
    let mut x = left.clone();

    for (a, b) in x.as_mut().iter_mut().zip(right.as_ref().iter()) {
        *a = a.wrapping_add(b.wrapping_mul(3u64));
    }

    let mut out = LweCiphertext::new(0u64, left.lwe_size(), ciphertext_modulus);

    pbs_first_nonzero(
        &x,
        ksk,
        bsk,
        lwe_secret_key_after_ks,
        polynomial_size,
        glwe_size,
        ciphertext_modulus,
        delta,
        &mut out,
    );

    out
}

#[allow(non_snake_case)]
pub fn sgn(
    T: &[LweCiphertext<Vec<u64>>],
    P: &[LweCiphertext<Vec<u64>>],
    parts_per_char: usize,
    ciphertext_modulus: CiphertextModulus<u64>,
    lwe_secret_key_after_ks: &LweSecretKey<Vec<u64>>,
    lwe_secret_key: &LweSecretKey<Vec<u64>>,
    ksk: &LweKeyswitchKey<Vec<u64>>,
    bsk: &FourierLweBootstrapKey<aligned_vec::ABox<[tfhe::core_crypto::fft_impl::fft64::c64]>>,
    polynomial_size: PolynomialSize,
    glwe_size: GlweSize,
    delta: u64,
) -> LweCiphertext<Vec<u64>> {
    assert_eq!(T.len(), P.len());
    assert!(!T.is_empty());
    assert!(parts_per_char > 0);
    assert_eq!(T.len() % parts_per_char, 0);

    let lwe_size = lwe_secret_key.lwe_dimension().to_lwe_size();
    let m = T.len() / parts_per_char;

    //
    // Phase 1:
    // Compute the comparison result for each character independently. A
    // character spans `parts_per_char` LWE ciphertexts (radix-4 parts,
    // ordered MSB -> LSB); they are compared part-by-part and then combined
    // lexicographically (leftmost, i.e. most-significant, non-zero part
    // wins) before the pattern-side wildcard mask is applied. For
    // parts_per_char == 1 this reduces to a single-part comparison, matching
    // the previous behavior exactly.
    //
    // cmp_i:
    //   0  : equal or wildcard
    //   1  : T[i] > P[i]
    //   15 : T[i] < P[i] (= -1 mod 16)
    //
    let mut cmp_vec: Vec<LweCiphertext<Vec<u64>>> = (0..m)
        .into_par_iter()
        .map(|char_idx| {
            let base = char_idx * parts_per_char;

            // Raw per-part comparison, without the wildcard mask applied yet.
            let mut raw_cmp_parts: Vec<LweCiphertext<Vec<u64>>> = (0..parts_per_char)
                .map(|k| {
                    let ct_t = &T[base + k];
                    let ct_p = &P[base + k];

                    // diff = T_part - P_part
                    let mut diff = ct_t.clone();

                    for (d, p) in diff.as_mut().iter_mut().zip(ct_p.as_ref().iter()) {
                        *d = d.wrapping_sub(*p);
                    }

                    // cmp_raw:
                    //   0  if equal
                    //   1  if T_part > P_part
                    //   15 if T_part < P_part
                    let mut cmp_raw = LweCiphertext::new(0u64, lwe_size, ciphertext_modulus);

                    pbs_cmp_pm1_zero(
                        &diff,
                        ksk,
                        bsk,
                        lwe_secret_key_after_ks,
                        polynomial_size,
                        glwe_size,
                        ciphertext_modulus,
                        delta,
                        &mut cmp_raw,
                    );

                    cmp_raw
                })
                .collect();

            // Combine parts MSB -> LSB: the leftmost non-zero part decides
            // the character-level comparison (lexicographic radix-4 compare).
            while raw_cmp_parts.len() > 1 {
                raw_cmp_parts = raw_cmp_parts
                    .chunks(2)
                    .map(|chunk| {
                        if chunk.len() == 1 {
                            chunk[0].clone()
                        } else {
                            combine_first_nonzero(
                                &chunk[0],
                                &chunk[1],
                                ciphertext_modulus,
                                lwe_secret_key_after_ks,
                                ksk,
                                bsk,
                                polynomial_size,
                                glwe_size,
                                delta,
                            )
                        }
                    })
                    .collect();
            }
            let combined_raw_cmp = raw_cmp_parts.pop().unwrap();

            // wild_bit = 1 iff the pattern character is a wildcard, i.e. its
            // most-significant part (part 0) is 0 -- by construction every
            // part of a padding character is 0.
            let mut wild_bit = LweCiphertext::new(0u64, lwe_size, ciphertext_modulus);

            pbs_is_zero(
                &P[base],
                ksk,
                bsk,
                lwe_secret_key_after_ks,
                polynomial_size,
                glwe_size,
                ciphertext_modulus,
                delta,
                &mut wild_bit,
            );

            // x = cmp_raw + 3 * wild_bit
            //
            // wild=0:
            //   cmp=0   -> x=0  -> 0
            //   cmp=1   -> x=1  -> 1
            //   cmp=15  -> x=15 -> 15
            //
            // wild=1:
            //   cmp=0   -> x=3  -> 0
            //   cmp=1   -> x=4  -> 0
            //   cmp=15  -> x=2  -> 0
            let mut x = combined_raw_cmp;

            for (a, w) in x.as_mut().iter_mut().zip(wild_bit.as_ref().iter()) {
                *a = a.wrapping_add(w.wrapping_mul(3u64));
            }

            let mut cmp_i = LweCiphertext::new(0u64, lwe_size, ciphertext_modulus);

            pbs_mask_wild_cmp(
                &x,
                ksk,
                bsk,
                lwe_secret_key_after_ks,
                polynomial_size,
                glwe_size,
                ciphertext_modulus,
                delta,
                &mut cmp_i,
            );

            cmp_i
        })
        .collect();

    //
    // Phase 2:
    // Reduce to the leftmost non-zero comparison result.
    //
    // combine(a, b) = if a == 0 { b } else { a }
    //
    // Pairwise reduction preserves left-to-right order while allowing
    // all pairs at the same tree level to be evaluated in parallel.
    //
    while cmp_vec.len() > 1 {
        cmp_vec = cmp_vec
            .par_chunks(2)
            .map(|chunk| {
                if chunk.len() == 1 {
                    chunk[0].clone()
                } else {
                    combine_first_nonzero(
                        &chunk[0],
                        &chunk[1],
                        ciphertext_modulus,
                        lwe_secret_key_after_ks,
                        ksk,
                        bsk,
                        polynomial_size,
                        glwe_size,
                        delta,
                    )
                }
            })
            .collect();
    }

    cmp_vec.pop().unwrap()
}
