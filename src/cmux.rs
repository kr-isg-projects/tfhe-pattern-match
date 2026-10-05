use aligned_vec::ABox;
use std::collections::HashMap;

use tfhe::core_crypto::entities::FourierLweBootstrapKey;
use tfhe::core_crypto::fft_impl::fft64::c64;
use tfhe::core_crypto::prelude::*;

use crate::automorphism::AutomorphKey;
use crate::lwe_to_glwe::convert_lwe_to_glwe_by_trace_with_preprocessing;
use crate::pbs;

pub fn cmux_lwe_vec<C: Container<Element = c64>>(
    dst: &mut [LweCiphertext<Vec<u64>>],
    candidate: &[LweCiphertext<Vec<u64>>],
    selector: &FourierGgswCiphertext<C>,
    auto_keys: &HashMap<usize, AutomorphKey>,
    ciphertext_modulus: CiphertextModulus<u64>,
) {
    dst.iter_mut()
        .zip(candidate.iter())
        .for_each(|(dst_part, cand_part)| {
            let mut glwe_dst = GlweCiphertext::new(
                0u64,
                selector.glwe_size(),
                selector.polynomial_size(),
                ciphertext_modulus,
            );
            let mut glwe_cand = GlweCiphertext::new(
                0u64,
                selector.glwe_size(),
                selector.polynomial_size(),
                ciphertext_modulus,
            );
            convert_lwe_to_glwe_by_trace_with_preprocessing(dst_part, &mut glwe_dst, auto_keys);
            convert_lwe_to_glwe_by_trace_with_preprocessing(cand_part, &mut glwe_cand, auto_keys);
            cmux_assign(&mut glwe_dst, &mut glwe_cand, selector);

            let mut out_lwe = LweCiphertext::new(0u64, dst_part.lwe_size(), ciphertext_modulus);
            extract_lwe_sample_from_glwe_ciphertext(&glwe_dst, &mut out_lwe, MonomialDegree(0));
            *dst_part = out_lwe;
        });
}

#[allow(clippy::too_many_arguments)]
pub fn cmux_lwe_vec_with_refresh<C: Container<Element = c64>>(
    dst: &mut [LweCiphertext<Vec<u64>>],
    candidate: &[LweCiphertext<Vec<u64>>],
    selector: &FourierGgswCiphertext<C>,
    auto_keys: &HashMap<usize, AutomorphKey>,
    ksk: &LweKeyswitchKey<Vec<u64>>,
    bsk: &FourierLweBootstrapKey<ABox<[c64]>>,
    lwe_secret_key_after_ks: &LweSecretKey<Vec<u64>>,
    polynomial_size: PolynomialSize,
    glwe_size: GlweSize,
    ciphertext_modulus: CiphertextModulus<u64>,
    delta: u64,
) {
    dst.iter_mut()
        .zip(candidate.iter())
        .for_each(|(dst_part, cand_part)| {
            let mut glwe_dst =
                GlweCiphertext::new(0u64, glwe_size, polynomial_size, ciphertext_modulus);
            let mut glwe_cand =
                GlweCiphertext::new(0u64, glwe_size, polynomial_size, ciphertext_modulus);
            convert_lwe_to_glwe_by_trace_with_preprocessing(dst_part, &mut glwe_dst, auto_keys);
            convert_lwe_to_glwe_by_trace_with_preprocessing(cand_part, &mut glwe_cand, auto_keys);
            cmux_assign(&mut glwe_dst, &mut glwe_cand, selector);

            let mut out_lwe = LweCiphertext::new(0u64, dst_part.lwe_size(), ciphertext_modulus);
            let mut refreshed = LweCiphertext::new(0u64, dst_part.lwe_size(), ciphertext_modulus);
            extract_lwe_sample_from_glwe_ciphertext(&glwe_dst, &mut out_lwe, MonomialDegree(0));
            pbs::pbs_identity(
                &out_lwe,
                ksk,
                bsk,
                lwe_secret_key_after_ks,
                polynomial_size,
                glwe_size,
                ciphertext_modulus,
                delta,
                &mut refreshed,
            );
            *dst_part = refreshed;
        });
}
