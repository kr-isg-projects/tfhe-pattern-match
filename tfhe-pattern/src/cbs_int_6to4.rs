use tfhe::core_crypto::prelude::*;
use tfhe::core_crypto::fft_impl::fft64::c64;
use tfhe::core_crypto::entities::FourierLweBootstrapKey;
use refined_tfhe_lhe::{
    AutomorphKey,
    blind_rotate_for_msb,
    glwe_ciphertext_clone_from,
    convert_to_ggsw_after_blind_rotate,
    glwe_ciphertext_monic_monomial_div,
};
use crate::param::HYBRID_BASE_64;
use aligned_vec::ABox;
use std::collections::HashMap;

pub fn cbs_int_6to4(
    lwe: &LweCiphertext<Vec<u64>>,
    ksk: &LweKeyswitchKey<Vec<u64>>,
    bsk: &FourierLweBootstrapKey<ABox<[c64]>>, 
    auto_keys: &HashMap<usize, AutomorphKey<ABox<[c64]>>>,
    ss_key: &FourierGgswCiphertextList<Vec<c64>>,
    num_ggsw: usize,
) -> GgswCiphertextList<Vec<u64>> {
    let param = *HYBRID_BASE_64;
    let extract_size = 1;
    let message_size = 6;
    let glwe_dimension = param.glwe_dimension();
    let polynomial_size = param.polynomial_size();
    let cbs_base_log = param.cbs_base_log();
    let cbs_level = param.cbs_level();
    let ciphertext_modulus = param.ciphertext_modulus();
    let glwe_size = glwe_dimension.to_glwe_size();
    let mut buf = LweCiphertext::new(u64::ZERO, lwe.lwe_size(), ciphertext_modulus);
    buf.as_mut().clone_from_slice(lwe.as_ref());
    let mut ggsw_list_out = GgswCiphertextList::new(
        0u64,
        glwe_size,
        polynomial_size,
        cbs_base_log,
        cbs_level,
        GgswCiphertextCount(num_ggsw),
        ciphertext_modulus,
    );
    for (idx, mut ggsw_chunk) in ggsw_list_out.chunks_exact_mut(extract_size).enumerate() {
        let mut acc_glev = GlweCiphertextList::new(
            u64::ZERO,
            glwe_size,
            polynomial_size,
            GlweCiphertextCount(cbs_level.0),
            ciphertext_modulus,
        );
        let mut lwe_extract = LweCiphertext::new(u64::ZERO, buf.lwe_size(), ciphertext_modulus);
        lwe_ciphertext_cleartext_mul(
            &mut lwe_extract,
            &buf,
            Cleartext(1u64 << (message_size - extract_size * (idx + 1))),
        );
        let mut lwe_extract_ks = LweCiphertext::new(u64::ZERO, ksk.output_lwe_size(), ciphertext_modulus);
        keyswitch_lwe_ciphertext(
            ksk,
            &lwe_extract,
            &mut lwe_extract_ks,
        );
        blind_rotate_for_msb(
            &lwe_extract_ks,
            &mut acc_glev,
            bsk.as_view(),
            param.log_lut_count(),
            cbs_base_log,
            cbs_level,
            extract_size,
            ciphertext_modulus,
        );
        let mut fourier_ggsw_chunk_out = FourierGgswCiphertextList::new(
            vec![c64::default();
                extract_size * polynomial_size.to_fourier_polynomial_size().0
                    * glwe_size.0
                    * glwe_size.0
                    * cbs_level.0
            ],
            extract_size,
            glwe_size,
            polynomial_size,
            cbs_base_log,
            cbs_level,
        );
        let mut buf_next = LweCiphertext::new(u64::ZERO, buf.lwe_size(), ciphertext_modulus);
        let log_scale = u64::BITS as usize - message_size + idx * extract_size;
        let acc_plaintext = PlaintextList::from_container((0..polynomial_size.0).map(|i| {
            if i < (1 << extract_size) {
                if (i >> (extract_size - 1)) == 0 {
                    (i << log_scale) as u64
                } else {
                    (((1 << (extract_size - 1)) + ((1 << extract_size) - 1 - i)) << log_scale) as u64
                }
            } else {
                u64::ZERO
            }
        }).collect::<Vec<u64>>());
        let acc_id = allocate_and_trivially_encrypt_new_glwe_ciphertext(
            glwe_size,
            &acc_plaintext,
            ciphertext_modulus,
        );
        let mut ct0 = GlweCiphertext::new(u64::ZERO, glwe_size, polynomial_size, ciphertext_modulus);
        let mut ct1 = GlweCiphertext::new(u64::ZERO, glwe_size, polynomial_size, ciphertext_modulus);
        glwe_ciphertext_clone_from(
            &mut ct0,
            &acc_id,
        );
        for (i, (mut ggsw, mut fourier_ggsw)) in ggsw_chunk.iter_mut()
        .zip(fourier_ggsw_chunk_out.as_mut_view().into_ggsw_iter())
        .enumerate()
        {
            convert_to_ggsw_after_blind_rotate(
                &acc_glev,
                &mut ggsw,
                extract_size - i - 1,
                auto_keys,
                ss_key.as_view(),
                ciphertext_modulus,
            );
            convert_standard_ggsw_ciphertext_to_fourier(
                &ggsw,
                &mut fourier_ggsw,
            );
            glwe_ciphertext_monic_monomial_div(
                &mut ct1,
                &ct0,
                MonomialDegree(1 << i),
            );
            cmux_assign(
                &mut ct0,
                &mut ct1,
                &fourier_ggsw,
            );
        }
        let mut lwe_extract = LweCiphertext::new(u64::ZERO, lwe.lwe_size(), ciphertext_modulus);
        extract_lwe_sample_from_glwe_ciphertext(
            &ct0,
            &mut lwe_extract,
            MonomialDegree(0),
        );
        lwe_ciphertext_sub(
            &mut buf_next,
            &buf,
            &lwe_extract,
        );
        buf.as_mut().clone_from_slice(buf_next.as_ref());
    }
    ggsw_list_out
}

pub fn cbs_int_6to1(
    lwe: &LweCiphertext<Vec<u64>>,
    ksk: &LweKeyswitchKey<Vec<u64>>,
    bsk: &FourierLweBootstrapKey<ABox<[c64]>>,
    auto_keys: &HashMap<usize, AutomorphKey<ABox<[c64]>>>,
    ss_key: &FourierGgswCiphertextList<Vec<c64>>,
    idx: usize,
) -> GgswCiphertext<Vec<u64>> {
    let param = *HYBRID_BASE_64;
    let extract_size = 1;
    let message_size = 6;
    let glwe_dimension = param.glwe_dimension();
    let polynomial_size = param.polynomial_size();
    let cbs_base_log = param.cbs_base_log();
    let cbs_level = param.cbs_level();
    let ciphertext_modulus = param.ciphertext_modulus();
    let glwe_size = glwe_dimension.to_glwe_size();
    let mut buf = LweCiphertext::new(u64::ZERO, lwe.lwe_size(), ciphertext_modulus);
    buf.as_mut().clone_from_slice(lwe.as_ref());
    let mut ggsw = GgswCiphertext::new(
        0u64,
        glwe_size,
        polynomial_size,
        cbs_base_log,
        cbs_level,
        ciphertext_modulus,
    );
    for i in 0..=idx {
        let mut acc_glev = GlweCiphertextList::new(
            u64::ZERO,
            glwe_size,
            polynomial_size,
            GlweCiphertextCount(cbs_level.0),
            ciphertext_modulus,
        );
        let mut lwe_extract = LweCiphertext::new(u64::ZERO, buf.lwe_size(), ciphertext_modulus);
        lwe_ciphertext_cleartext_mul(
            &mut lwe_extract,
            &buf,
            Cleartext(1u64 << (message_size - extract_size * (i + 1))),
        );
        let mut lwe_extract_ks = LweCiphertext::new(u64::ZERO, ksk.output_lwe_size(), ciphertext_modulus);
        keyswitch_lwe_ciphertext(
            ksk,
            &lwe_extract,
            &mut lwe_extract_ks,
        );
        blind_rotate_for_msb(
            &lwe_extract_ks,
            &mut acc_glev,
            bsk.as_view(),
            param.log_lut_count(),
            cbs_base_log,
            cbs_level,
            extract_size,
            ciphertext_modulus,
        );
        let log_scale = u64::BITS as usize - message_size + i * extract_size;
        let acc_plaintext = PlaintextList::from_container((0..polynomial_size.0).map(|j| {
            if j < (1 << extract_size) {
                if (j >> (extract_size - 1)) == 0 {
                    (j << log_scale) as u64
                } else {
                    (((1 << (extract_size - 1)) + ((1 << extract_size) - 1 - j)) << log_scale) as u64
                }
            } else {
                u64::ZERO
            }
        }).collect::<Vec<u64>>());
        let acc_id = allocate_and_trivially_encrypt_new_glwe_ciphertext(
            glwe_size,
            &acc_plaintext,
            ciphertext_modulus,
        );
        let mut ct0 = GlweCiphertext::new(u64::ZERO, glwe_size, polynomial_size, ciphertext_modulus);
        glwe_ciphertext_clone_from(&mut ct0, &acc_id);
        convert_to_ggsw_after_blind_rotate(
            &acc_glev,
            &mut ggsw,
            0,
            auto_keys,
            ss_key.as_view(),
            ciphertext_modulus,
        );
        if i == idx {
            break;
        }
        let mut lwe_extract2 = LweCiphertext::new(u64::ZERO, buf.lwe_size(), ciphertext_modulus);
        extract_lwe_sample_from_glwe_ciphertext(
            &ct0,
            &mut lwe_extract2,
            MonomialDegree(0),
        );
        let mut buf_next = LweCiphertext::new(u64::ZERO, buf.lwe_size(), ciphertext_modulus);
        lwe_ciphertext_sub(&mut buf_next, &buf, &lwe_extract2);
        buf.as_mut().clone_from_slice(buf_next.as_ref());
    }
    ggsw
}