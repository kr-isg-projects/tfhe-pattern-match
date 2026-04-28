use tfhe::core_crypto::entities::{LweCiphertext, LweSecretKey};
use tfhe::core_crypto::prelude::*;
use rayon::prelude::*;
use crate::pbs::*;

#[allow(non_snake_case)]
pub fn sgn(
    T: &[LweCiphertext<Vec<u64>>],
    P: &[LweCiphertext<Vec<u64>>],
    ciphertext_modulus: CiphertextModulus<u64>,
    lwe_secret_key_after_ks: &LweSecretKey<Vec<u64>>,
    lwe_secret_key: &LweSecretKey<Vec<u64>>,
    ksk: &LweKeyswitchKey<Vec<u64>>,
    bsk: &tfhe::core_crypto::entities::FourierLweBootstrapKey<aligned_vec::ABox<[tfhe::core_crypto::fft_impl::fft64::c64]>>,
    polynomial_size2: PolynomialSize,
    glwe_size: GlweSize,
    delta: u64,
) -> LweCiphertext<Vec<u64>> {
    let mut sgn_vec: Vec<LweCiphertext<Vec<u64>>> = T.par_iter().zip(P.par_iter())
        .map(|(ct0, ct1)| {
            let mut diff = ct0.clone();
            for (d, c1) in diff.as_mut().iter_mut().zip(ct1.as_ref().iter()) {
                *d = d.wrapping_sub(*c1);
            }
            let mut out = LweCiphertext::new(0u64, lwe_secret_key.lwe_dimension().to_lwe_size(), ciphertext_modulus);
            pbs_sgn_with_special_words(
                &diff,
                ksk,
                bsk,
                lwe_secret_key_after_ks,
                polynomial_size2,
                glwe_size,
                ciphertext_modulus,
                delta,
                &mut out,
            );
            out
        })
        .collect();

    while sgn_vec.len() > 1 {
        sgn_vec = sgn_vec.par_chunks(5).map(|chunk| {
            let a = &chunk[0];
            let b = if chunk.len() > 1 { &chunk[1] } else { &chunk[0] };
            let c = if chunk.len() > 2 { &chunk[2] } else { &chunk[0] };
            let d = if chunk.len() > 3 { &chunk[3] } else { &chunk[0] };
            let e = if chunk.len() > 4 { &chunk[4] } else { &chunk[0] };

            let mut val = a.clone();
            for v in val.as_mut().iter_mut() {
                *v = v.wrapping_mul(16);
            }
            let mut tmp = b.clone();
            for v in tmp.as_mut().iter_mut() {
                *v = v.wrapping_mul(8);
            }
            for (v, t) in val.as_mut().iter_mut().zip(tmp.as_ref().iter()) {
                *v = v.wrapping_add(*t);
            }
            let mut tmp = c.clone();
            for v in tmp.as_mut().iter_mut() {
                *v = v.wrapping_mul(4);
            }
            for (v, t) in val.as_mut().iter_mut().zip(tmp.as_ref().iter()) {
                *v = v.wrapping_add(*t);
            }
            let mut tmp = d.clone();
            for v in tmp.as_mut().iter_mut() {
                *v = v.wrapping_mul(2);
            }
            for (v, t) in val.as_mut().iter_mut().zip(tmp.as_ref().iter()) {
                *v = v.wrapping_add(*t);
            }
            for (v, t) in val.as_mut().iter_mut().zip(e.as_ref().iter()) {
                *v = v.wrapping_add(*t);
            }

            let mut out = LweCiphertext::new(0u64, val.lwe_size(), ciphertext_modulus);
            pbs_sgn(
                &val,
                ksk,
                bsk,
                lwe_secret_key_after_ks,
                polynomial_size2,
                glwe_size,
                ciphertext_modulus,
                delta,
                &mut out,
            );
            out
        }).collect();
    }
    sgn_vec[0].clone()
}