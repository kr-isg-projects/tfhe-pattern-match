pub mod arithmetic;
pub mod automorphism;
pub mod cbs;
pub mod cmux;
pub mod keygen;
pub mod lut_eval;
pub mod lwe_to_glwe;
pub mod nfbe;
pub mod param;
pub mod pbs;
pub mod sgn;
pub mod utils;
pub mod wop_pbs;

use crate::arithmetic::ArithmeticCtx;
use crate::cmux as cmux_helpers;
use crate::param::{FheKeys, FheParams};
use rayon::prelude::*;
use tfhe::core_crypto::prelude::*;

pub struct FheFindResult {
    pub found: LweCiphertext<Vec<u64>>,
    pub index: Vec<LweCiphertext<Vec<u64>>>,
}

pub fn fhe_find(
    ct_text: Vec<GlweCiphertextList<Vec<u64>>>,
    ct_text_shifted: Vec<GlweCiphertextList<Vec<u64>>>,
    ct_sa: Vec<GlweCiphertextList<Vec<u64>>>,
    ct_pattern: Vec<LweCiphertext<Vec<u64>>>,
    keys: FheKeys,
    params: FheParams,
) -> FheFindResult {
    let parts_per_char = ct_text.len();
    assert_eq!(ct_text_shifted.len(), parts_per_char);
    assert_eq!(ct_pattern.len(), params.m * parts_per_char);
    let FheKeys {
        lwe_secret_key,
        glwe_secret_key,
        lwe_secret_key_after_ks,
        bsk,
        ksk,
        pfpksk,
        auto_keys,
    } = keys;
    let FheParams {
        log_n,
        n,
        m,
        n_parts,
        top_part_bits,
        polynomial_size,
        glwe_size,
        glwe_modular_std_dev,
        cbs_base_log,
        cbs_level,
        ciphertext_modulus,
        mut encryption_generator,
        delta,
    } = params;
    let noise_distribution = Gaussian::from_dispersion_parameter(glwe_modular_std_dev, 0.0);
    let l_val = 0u64;
    let r_val = n - 1;
    let c_val = (n - 1) / 2;
    let tmp = LweCiphertext::new(
        0u64,
        lwe_secret_key.lwe_dimension().to_lwe_size(),
        ciphertext_modulus,
    );
    let mut i_parts: Vec<LweCiphertext<Vec<u64>>> = {
        let mut v = Vec::with_capacity(n_parts);
        for _ in 0..n_parts {
            v.push(tmp.clone());
        }
        v
    };
    let c_val_bits: Vec<usize> = (0..log_n).map(|i| ((c_val >> i) & 1) as usize).collect();
    let mut ggsw_list_c = GgswCiphertextList::new(
        0u64,
        glwe_size,
        polynomial_size,
        cbs_base_log,
        cbs_level,
        GgswCiphertextCount(log_n),
        ciphertext_modulus,
    );
    for (bit, mut ggsw) in c_val_bits.iter().rev().zip(ggsw_list_c.iter_mut()) {
        encrypt_constant_ggsw_ciphertext(
            &glwe_secret_key,
            &mut ggsw,
            Cleartext(*bit as u64),
            noise_distribution,
            &mut encryption_generator,
        );
    }
    let mut l_parts = utils::init_constant_parts_lwe(
        l_val,
        n_parts,
        &lwe_secret_key,
        ciphertext_modulus,
        noise_distribution,
        &mut encryption_generator,
        delta,
    );
    let mut r_parts = utils::init_radix_parts_lwe(
        r_val,
        n_parts,
        &lwe_secret_key,
        ciphertext_modulus,
        noise_distribution,
        &mut encryption_generator,
        delta,
    );
    let mut c_parts = utils::init_radix_parts_lwe(
        c_val,
        n_parts,
        &lwe_secret_key,
        ciphertext_modulus,
        noise_distribution,
        &mut encryption_generator,
        delta,
    );
    let mut found = LweCiphertext::new(
        0u64,
        lwe_secret_key.lwe_dimension().to_lwe_size(),
        ciphertext_modulus,
    );
    let profile_enabled = std::env::var_os("FHE_FIND_PROFILE").is_some();
    let arithmetic_ctx = ArithmeticCtx {
        lwe_secret_key: &lwe_secret_key,
        lwe_secret_key_after_ks: &lwe_secret_key_after_ks,
        ksk: &ksk,
        bsk: &bsk,
        polynomial_size,
        glwe_size,
        ciphertext_modulus,
        delta,
    };
    for iter in 0..log_n {
        let part_build_start = std::time::Instant::now();
        for i in 0..n_parts {
            let lwe_out = lut_eval::private_lut_eval_n_to_4(
                &ct_sa[i],
                &ggsw_list_c,
                &lwe_secret_key,
                polynomial_size,
                cbs_base_log,
                cbs_level,
                ciphertext_modulus,
            );
            i_parts[i] = lwe_out;
        }
        if profile_enabled {
            eprintln!(
                "[fhe_find iter={iter}] part_build: {:.3}s",
                part_build_start.elapsed().as_secs_f64(),
            );
        }

        let true_slice_start = std::time::Instant::now();
        let ct_true_slice = nfbe::extract_subsequence(
            &ct_text,
            &ct_text_shifted,
            &i_parts,
            m,
            n_parts,
            top_part_bits,
            &pfpksk,
            cbs_base_log,
            cbs_level,
            &arithmetic_ctx,
        );
        if profile_enabled {
            eprintln!(
                "[fhe_find iter={iter}] true_slice: {:.3}s",
                true_slice_start.elapsed().as_secs_f64(),
            );
        }

        let sgn_start = std::time::Instant::now();
        let mut actual_sgn = sgn::sgn(
            &ct_true_slice,
            &ct_pattern,
            parts_per_char,
            ciphertext_modulus,
            &lwe_secret_key_after_ks,
            &lwe_secret_key,
            &ksk,
            &bsk,
            polynomial_size,
            glwe_size,
            delta,
        );
        if profile_enabled {
            eprintln!(
                "[fhe_find iter={iter}] sgn: {:.3}s",
                sgn_start.elapsed().as_secs_f64(),
            );
        }
        let post_sgn_cbs_start = std::time::Instant::now();
        let lt = cbs::cbs_int_4to1(&actual_sgn, &ksk, &bsk, &pfpksk, 1);
        let lwe_plain1 = utils::encrypt_small_int(
            1,
            &lwe_secret_key,
            ciphertext_modulus,
            noise_distribution,
            &mut encryption_generator,
            delta,
        );
        for (x, y) in actual_sgn
            .as_mut()
            .iter_mut()
            .zip(lwe_plain1.as_ref().iter())
        {
            *x = x.wrapping_add(*y);
        }
        if iter == log_n - 1 {
            pbs::pbs_eq_one(
                &actual_sgn,
                &ksk,
                &bsk,
                &lwe_secret_key_after_ks,
                polynomial_size,
                glwe_size,
                ciphertext_modulus,
                delta,
                &mut found,
            );
            break;
        }
        let gt = cbs::cbs_int_4to1(&actual_sgn, &ksk, &bsk, &pfpksk, 1);
        if profile_enabled {
            eprintln!(
                "[fhe_find iter={iter}] post_sgn_cbs: {:.3}s",
                post_sgn_cbs_start.elapsed().as_secs_f64(),
            );
        }

        let c_pm1_start = std::time::Instant::now();
        let c_parts_p1 = arithmetic::add_one_lwe_vec(
            &c_parts,
            noise_distribution,
            &mut encryption_generator,
            &arithmetic_ctx,
        );
        let c_parts_m1 = arithmetic::sub_one_lwe_vec(
            &c_parts,
            noise_distribution,
            &mut encryption_generator,
            &arithmetic_ctx,
        );
        if profile_enabled {
            eprintln!(
                "[fhe_find iter={iter}] c_pm1: {:.3}s",
                c_pm1_start.elapsed().as_secs_f64(),
            );
        }

        let lr_update_start = std::time::Instant::now();
        let fourier_ggsw_lt = {
            let mut fourier =
                FourierGgswCiphertext::new(glwe_size, polynomial_size, cbs_base_log, cbs_level);
            convert_standard_ggsw_ciphertext_to_fourier(&lt, &mut fourier);
            fourier
        };
        cmux_helpers::cmux_lwe_vec(
            &mut l_parts,
            &c_parts_p1,
            &fourier_ggsw_lt,
            &auto_keys,
            ciphertext_modulus,
        );
        let fourier_ggsw_gt = {
            let mut fourier =
                FourierGgswCiphertext::new(glwe_size, polynomial_size, cbs_base_log, cbs_level);
            convert_standard_ggsw_ciphertext_to_fourier(&gt, &mut fourier);
            fourier
        };
        cmux_helpers::cmux_lwe_vec(
            &mut r_parts,
            &c_parts_m1,
            &fourier_ggsw_gt,
            &auto_keys,
            ciphertext_modulus,
        );
        if profile_enabled {
            eprintln!(
                "[fhe_find iter={iter}] lr_update: {:.3}s",
                lr_update_start.elapsed().as_secs_f64(),
            );
        }

        let c_sum_start = std::time::Instant::now();
        let sum_parts = arithmetic::add_lwe_vec(
            &l_parts,
            &r_parts,
            noise_distribution,
            &mut encryption_generator,
            &arithmetic_ctx,
        );
        if profile_enabled {
            eprintln!(
                "[fhe_find iter={iter}] c_sum: {:.3}s",
                c_sum_start.elapsed().as_secs_f64(),
            );
        }

        let c_div2_start = std::time::Instant::now();
        c_parts = arithmetic::div2_lwe_vec(
            &sum_parts,
            n_parts,
            noise_distribution,
            &mut encryption_generator,
            &arithmetic_ctx,
        );
        if profile_enabled {
            eprintln!(
                "[fhe_find iter={iter}] c_div2: {:.3}s",
                c_div2_start.elapsed().as_secs_f64(),
            );
        }

        let c_cbs_start = std::time::Instant::now();
        let ggsw_lists_c: Vec<_> = c_parts
            .par_iter()
            .enumerate()
            .map(|(i, ct)| {
                let bits = if i == n_parts - 1 { top_part_bits } else { 2 };

                cbs::cbs_int_4to2(ct, &ksk, &bsk, &pfpksk, bits)
            })
            .collect();
        ggsw_list_c = utils::concat_ggsw_lists(&ggsw_lists_c);
        if profile_enabled {
            eprintln!(
                "[fhe_find iter={iter}] c_cbs: {:.3}s",
                c_cbs_start.elapsed().as_secs_f64(),
            );
        }
    }
    FheFindResult {
        found,
        index: i_parts,
    }
}
