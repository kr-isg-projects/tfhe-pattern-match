
pub mod cbs_int_6to4;
pub mod param;
pub mod lut_eval;
pub mod pbs;
pub mod sgn;

use tfhe::core_crypto::prelude::*;
use tfhe::core_crypto::algorithms::encrypt_lwe_ciphertext;
use tfhe::core_crypto::entities::GlweCiphertext;
use refined_tfhe_lhe::convert_lwe_to_glwe_by_trace_with_preprocessing;
use crate::param::{FheKeys, FheParams};

fn concat_ggsw_lists(
    lists: &[GgswCiphertextList<Vec<u64>>]
) -> GgswCiphertextList<Vec<u64>> {
    let mut all_ggsw = Vec::new();
    for ggsw_list in lists.iter().rev() {
        for ggsw in ggsw_list.iter().rev() {
            all_ggsw.push(ggsw.clone().into_container());
        }
    }
    let first = &lists[0];
    GgswCiphertextList::from_container(
        all_ggsw.concat(),
        first.glwe_size(),
        first.polynomial_size(),
        first.decomposition_base_log(),
        first.decomposition_level_count(),
        first.ciphertext_modulus(),
    )
}

pub fn fhe_find(
    ct_text: GlweCiphertextList<Vec<u64>>,
    ct_text_shifted: GlweCiphertextList<Vec<u64>>,
    ct_sa: Vec<GlweCiphertextList<Vec<u64>>>,
    ct_pattern: Vec<LweCiphertext<Vec<u64>>>,
    keys: FheKeys,
    params: FheParams,
) -> Vec<LweCiphertext<Vec<u64>>> {
    let FheKeys {
        lwe_secret_key,
        glwe_secret_key,
        lwe_secret_key_after_ks,
        bsk,
        ksk,
        auto_keys,
        ss_key,
    } = keys;
    let FheParams {
        log_n,
        n,
        m,
        n_nibbles,
        n_rems,
        polynomial_size,
        glwe_size,
        glwe_modular_std_dev,
        cbs_base_log,
        cbs_level,
        ciphertext_modulus,
        mut encryption_generator,
        delta,
    } = params;
    let l_val = 0u64;
    let r_val = n - 1;
    let c_val = (n - 1) / 2;
    let tmp = LweCiphertext::new(
        0u64,
        lwe_secret_key.lwe_dimension().to_lwe_size(),
        ciphertext_modulus,
    );
    let mut lwe_out_carry = tmp.clone();
    let mut i_parts: Vec<LweCiphertext<Vec<u64>>> = {
        let mut v = Vec::with_capacity(n_nibbles);
        for _ in 0..n_nibbles{
            v.push(tmp.clone());
        }
        v
    };
    let c_val_bits: Vec<usize> = (0..log_n)
        .map(|i| ((c_val >> i) & 1) as usize)
        .collect();
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
            Plaintext(*bit as u64),
            glwe_modular_std_dev,
            &mut encryption_generator,
        );
    }
    let mut l_parts = vec![];
    let mut r_parts = vec![];
    let mut c_parts = vec![];
    for _ in 0..n_nibbles{
        let mut ct = LweCiphertext::new(
            0u64,
            lwe_secret_key.lwe_dimension().to_lwe_size(),
            ciphertext_modulus,
        );
        encrypt_lwe_ciphertext(
            &lwe_secret_key,
            &mut ct,
            Plaintext(l_val * delta),
            glwe_modular_std_dev,
            &mut encryption_generator,
        );
        l_parts.push(ct);
    }
    for i in 0..n_nibbles{
        let part_val = (r_val >> (i * 4)) & 0xF;
        let mut ct = LweCiphertext::new(
            0u64,
            lwe_secret_key.lwe_dimension().to_lwe_size(),
            ciphertext_modulus,
        );
        encrypt_lwe_ciphertext(
            &lwe_secret_key,
            &mut ct,
            Plaintext(part_val * delta),
            glwe_modular_std_dev,
            &mut encryption_generator,
        );
        r_parts.push(ct);
    }
    for i in 0..n_nibbles{
        let part_val = (c_val >> (i * 4)) & 0xF;
        let mut ct = LweCiphertext::new(
            0u64,
            lwe_secret_key.lwe_dimension().to_lwe_size(),
            ciphertext_modulus,
        );
        encrypt_lwe_ciphertext(
            &lwe_secret_key,
            &mut ct,
            Plaintext(part_val * delta),
            glwe_modular_std_dev,
            &mut encryption_generator,
        );
        c_parts.push(ct);
    }
    let mut ggsw_lists_i: Vec<GgswCiphertextList<Vec<u64>>> = Vec::with_capacity(log_n);
    let mut ggsw_lists_i_shifted: Vec<GgswCiphertextList<Vec<u64>>> = Vec::with_capacity(log_n);
    let mut eq = GgswCiphertext::new(
        0u64,
        glwe_size,
        polynomial_size,
        cbs_base_log,
        cbs_level,
        ciphertext_modulus,
    );
    for iter in 0..log_n {
        ggsw_lists_i.clear();
        ggsw_lists_i_shifted.clear();
        for i in 0..n_nibbles{
            let lwe_out = lut_eval::private_lut_eval_n_to_4(
                &ct_sa[i],
                &ggsw_list_c,
                &lwe_secret_key,
                polynomial_size,
                cbs_base_log,
                cbs_level,
                ciphertext_modulus,
            );
            i_parts[i] = lwe_out.clone();
            let mut lwe_out_shifted = lwe_out.clone();
            if i == 2 {
                let lwe_plain8 = {
                    let mut ct = LweCiphertext::new(
                        0u64,
                        lwe_secret_key.lwe_dimension().to_lwe_size(),
                        ciphertext_modulus,
                    );
                    encrypt_lwe_ciphertext(
                        &lwe_secret_key,
                        &mut ct,
                        Plaintext(8 * delta),
                        glwe_modular_std_dev,
                        &mut encryption_generator,
                    );
                    ct
                };
                for (x, y) in lwe_out_shifted.as_mut().iter_mut().zip(lwe_plain8.as_ref().iter()) {
                    *x = x.wrapping_add(*y);
                }
                pbs::pbs_carry(
                    &lwe_out_shifted,
                    &ksk,
                    &bsk,
                    &lwe_secret_key_after_ks,
                    polynomial_size,
                    glwe_size,
                    ciphertext_modulus,
                    delta,
                    &mut lwe_out_carry,
                );
            } if i >= 3 {
                for (x, y) in lwe_out_shifted.as_mut().iter_mut().zip(lwe_out_carry.as_ref().iter()) {
                    *x = x.wrapping_add(*y);
                }
                pbs::pbs_carry(
                    &lwe_out_shifted,
                    &ksk,
                    &bsk,
                    &lwe_secret_key_after_ks,
                    polynomial_size,
                    glwe_size,
                    ciphertext_modulus,
                    delta,
                    &mut lwe_out_carry,
                );
            }
            let bits = if i == n_nibbles - 1 { n_rems } else { 4 };
            let ggsw_list_out = cbs_int_6to4::cbs_int_6to4(
                &lwe_out,
                &ksk,
                &bsk,
                &auto_keys,
                &ss_key,
                bits,
            );
            let ggsw_list_out_shifted = cbs_int_6to4::cbs_int_6to4(
                &lwe_out_shifted,
                &ksk,
                &bsk,
                &auto_keys,
                &ss_key,
                bits,
            );
            ggsw_lists_i.push(ggsw_list_out);
            ggsw_lists_i_shifted.push(ggsw_list_out_shifted);
        }
        let ggsw_list_i = concat_ggsw_lists(&ggsw_lists_i);
        let ggsw_list_i_shifted = concat_ggsw_lists(&ggsw_lists_i_shifted);

        let extract_len = (m as usize).min(polynomial_size.0 / 2);
        let ct_text_slice = lut_eval::private_lut_eval_n_to_4_slice_multi(
            &ct_text,
            &ggsw_list_i,
            &lwe_secret_key,
            polynomial_size,
            cbs_base_log,
            cbs_level,
            ciphertext_modulus,
            extract_len,
        );
        let ct_text_slice_shifted = lut_eval::private_lut_eval_n_to_4_slice_multi(
            &ct_text_shifted,
            &ggsw_list_i_shifted,
            &lwe_secret_key,
            polynomial_size,
            cbs_base_log,
            cbs_level,
            ciphertext_modulus,
            extract_len,
        );

        let lwe_sgn = sgn::sgn(
            &ct_text_slice,
            &ct_pattern,
            ciphertext_modulus,
            &lwe_secret_key_after_ks,
            &lwe_secret_key,
            &ksk,
            &bsk,
            polynomial_size,
            glwe_size,
            delta,
        );
        let lwe_sgn_shifted = sgn::sgn(
            &ct_text_slice_shifted,
            &ct_pattern,
            ciphertext_modulus,
            &lwe_secret_key_after_ks,
            &lwe_secret_key,
            &ksk,
            &bsk,
            polynomial_size,
            glwe_size,
            delta,
        );
        let mut lwe_sum = lwe_sgn.clone();
        for (s, b) in lwe_sum.as_mut().iter_mut().zip(lwe_sgn_shifted.as_ref().iter()) {
            let triple_b = b.wrapping_mul(3);
            *s = s.wrapping_add(triple_b);
        }
        let mut actual_sgn = LweCiphertext::new(0u64, lwe_sum.lwe_size(), ciphertext_modulus);
        pbs::pbs_true_sgn(
            &lwe_sum,
            &ksk,
            &bsk,
            &lwe_secret_key_after_ks,
            polynomial_size,
            glwe_size,
            ciphertext_modulus,
            delta,
            &mut actual_sgn,
        );

        let lt = cbs_int_6to4::cbs_int_6to1(
            &actual_sgn,
            &ksk,
            &bsk,
            &auto_keys,
            &ss_key,
            1,
        );
        let lwe_plain1 = {
            let mut ct = LweCiphertext::new(
                0u64,
                lwe_secret_key.lwe_dimension().to_lwe_size(),
                ciphertext_modulus,
            );
            encrypt_lwe_ciphertext(
                &lwe_secret_key,
                &mut ct,
                Plaintext(1 * delta),
                glwe_modular_std_dev,
                &mut encryption_generator,
            );
            ct
        };
        for (x, y) in actual_sgn.as_mut().iter_mut().zip(lwe_plain1.as_ref().iter()) {
            *x = x.wrapping_add(*y);
        }
        eq = cbs_int_6to4::cbs_int_6to1(
            &actual_sgn,
            &ksk,
            &bsk,
            &auto_keys,
            &ss_key,
            0,
        );
        if iter == log_n - 1 {
            break;
        }
        let gt = cbs_int_6to4::cbs_int_6to1(
            &actual_sgn,
            &ksk,
            &bsk,
            &auto_keys,
            &ss_key,
            1,
        );

        let mut c_parts_p1: Vec<LweCiphertext<Vec<u64>>> = Vec::with_capacity(c_parts.len());
        let mut c_parts_m1: Vec<LweCiphertext<Vec<u64>>> = Vec::with_capacity(c_parts.len());
        for ct in c_parts.iter() {
            let mut ct_copy = LweCiphertext::new(
                0u64,
                ct.lwe_size(), 
                ct.ciphertext_modulus(),
            );
            ct_copy.as_mut().copy_from_slice(ct.as_ref());
            c_parts_p1.push(ct_copy);
        }
        for ct in c_parts.iter() {
            let mut ct_copy = LweCiphertext::new(
                0u64,
                ct.lwe_size(), 
                ct.ciphertext_modulus(), 
            );
            ct_copy.as_mut().copy_from_slice(ct.as_ref());
            c_parts_m1.push(ct_copy);
        }

        let mut ct_carry_p = LweCiphertext::new(
            0u64,
            lwe_secret_key.lwe_dimension().to_lwe_size(),
            ciphertext_modulus,
        );
        let mut ct_tmp_p = LweCiphertext::new(
            0u64,
            lwe_secret_key.lwe_dimension().to_lwe_size(),
            ciphertext_modulus,
        );
        for (i, ct) in c_parts_p1.iter_mut().enumerate() {
            if i == 0 {
                encrypt_lwe_ciphertext(
                    &lwe_secret_key,
                    &mut ct_carry_p,
                    Plaintext(1 * delta),
                    glwe_modular_std_dev,
                    &mut encryption_generator,
                );
            }
            let x_slice = ct.as_ref();
            let y_slice = ct_carry_p.as_ref();
            let tmp_slice = ct_tmp_p.as_mut();
            tmp_slice.copy_from_slice(x_slice);
            for idx in 0..tmp_slice.len() {
                tmp_slice[idx] = tmp_slice[idx].wrapping_add(y_slice[idx]);
            }
            pbs::pbs_carry(
                &ct_tmp_p, 
                &ksk,
                &bsk,
                &lwe_secret_key_after_ks,
                polynomial_size,
                glwe_size,
                ciphertext_modulus,
                delta,
                &mut ct_carry_p,
            );
            pbs::pbs_reminder(
                &ct_tmp_p, 
                &ksk,
                &bsk,
                &lwe_secret_key_after_ks,
                polynomial_size,
                glwe_size,
                ciphertext_modulus,
                delta,
                ct, 
            );
        }
        let mut ct_carry = LweCiphertext::new(
            0u64,
            lwe_secret_key.lwe_dimension().to_lwe_size(),
            ciphertext_modulus,
        );
        let mut ct_tmp = LweCiphertext::new(
            0u64,
            lwe_secret_key.lwe_dimension().to_lwe_size(),
            ciphertext_modulus,
        );
        for (i, ct) in c_parts_m1.iter_mut().enumerate() {
            if i == 0 {
                encrypt_lwe_ciphertext(
                    &lwe_secret_key,
                    &mut ct_carry,
                    Plaintext(1 * delta),
                    glwe_modular_std_dev,
                    &mut encryption_generator,
                );
                let x_slice = ct.as_ref();
                let y_slice = ct_carry.as_ref();
                let tmp_slice = ct_tmp.as_mut();
                tmp_slice.copy_from_slice(x_slice);
                for idx in 0..tmp_slice.len() {
                    tmp_slice[idx] = tmp_slice[idx].wrapping_sub(y_slice[idx]);
                }
                pbs::pbs_carry_sub(
                    &ct_tmp, 
                    &ksk,
                    &bsk,
                    &lwe_secret_key_after_ks,
                    polynomial_size,
                    glwe_size,
                    ciphertext_modulus,
                    delta,
                    &mut ct_carry, 
                );
                pbs::pbs_reminder(
                    &ct_tmp, 
                    &ksk,
                    &bsk,
                    &lwe_secret_key_after_ks,
                    polynomial_size,
                    glwe_size,
                    ciphertext_modulus,
                    delta,
                    ct, 
                );
            } else {
                let x_slice = ct.as_ref();
                let y_slice = ct_carry.as_ref();
                let tmp_slice = ct_tmp.as_mut();
                tmp_slice.copy_from_slice(x_slice);
                for idx in 0..tmp_slice.len() {
                    tmp_slice[idx] = tmp_slice[idx].wrapping_sub(y_slice[idx]);
                }
                pbs::pbs_carry_sub(
                    &ct_tmp,
                    &ksk,
                    &bsk,
                    &lwe_secret_key_after_ks,
                    polynomial_size,
                    glwe_size,
                    ciphertext_modulus,
                    delta,
                    &mut ct_carry, 
                );
                pbs::pbs_reminder(
                    &ct_tmp,
                    &ksk,
                    &bsk,
                    &lwe_secret_key_after_ks,
                    polynomial_size,
                    glwe_size,
                    ciphertext_modulus,
                    delta,
                    ct,
                );
            }
        }
        let fourier_ggsw_lt = {
            let mut fourier = FourierGgswCiphertext::new(
                glwe_size,
                polynomial_size,
                cbs_base_log,
                cbs_level,
            );
            convert_standard_ggsw_ciphertext_to_fourier(&lt, &mut fourier);
            fourier
        };
        for i in 0..n_nibbles{
            let mut glwe_l = GlweCiphertext::new(0u64, glwe_size, polynomial_size, ciphertext_modulus);
            let mut glwe_cp1 = GlweCiphertext::new(0u64, glwe_size, polynomial_size, ciphertext_modulus);
            convert_lwe_to_glwe_by_trace_with_preprocessing(&l_parts[i], &mut glwe_l, &auto_keys);
            convert_lwe_to_glwe_by_trace_with_preprocessing(&c_parts_p1[i], &mut glwe_cp1, &auto_keys);
            cmux_assign(&mut glwe_l, &mut glwe_cp1, &fourier_ggsw_lt);
            let mut out_lwe = LweCiphertext::new(
                0u64,
                lwe_secret_key.lwe_dimension().to_lwe_size(),
                ciphertext_modulus,
            );
            let mut lwe_l = LweCiphertext::new(
                0u64,
                lwe_secret_key.lwe_dimension().to_lwe_size(),
                ciphertext_modulus,
            );
            extract_lwe_sample_from_glwe_ciphertext(&glwe_l, &mut out_lwe, MonomialDegree(0));
            pbs::pbs_identity(
                    &out_lwe,
                    &ksk,
                    &bsk,
                    &lwe_secret_key_after_ks,
                    polynomial_size,
                    glwe_size,
                    ciphertext_modulus,
                    delta,
                    &mut lwe_l,
                );
            l_parts[i].as_mut().copy_from_slice(lwe_l.as_ref());
        }
        let fourier_ggsw_gt = {
            let mut fourier = FourierGgswCiphertext::new(
                glwe_size,
                polynomial_size,
                cbs_base_log,
                cbs_level,
            );
            convert_standard_ggsw_ciphertext_to_fourier(&gt, &mut fourier);
            fourier
        };
        for i in 0..n_nibbles{
            let mut glwe_r = GlweCiphertext::new(0u64, glwe_size, polynomial_size, ciphertext_modulus);
            let mut glwe_cm1 = GlweCiphertext::new(0u64, glwe_size, polynomial_size, ciphertext_modulus);
            convert_lwe_to_glwe_by_trace_with_preprocessing(&r_parts[i], &mut glwe_r, &auto_keys);
            convert_lwe_to_glwe_by_trace_with_preprocessing(&c_parts_m1[i], &mut glwe_cm1, &auto_keys);
            cmux_assign(&mut glwe_r, &mut glwe_cm1, &fourier_ggsw_gt);
            let mut out_lwe = LweCiphertext::new(
                0u64,
                lwe_secret_key.lwe_dimension().to_lwe_size(),
                ciphertext_modulus,
            );
            let mut lwe_r = LweCiphertext::new(
                0u64,
                lwe_secret_key.lwe_dimension().to_lwe_size(),
                ciphertext_modulus,
            );
            extract_lwe_sample_from_glwe_ciphertext(&glwe_r, &mut out_lwe, MonomialDegree(0));
            pbs::pbs_identity(
                    &out_lwe,
                    &ksk,
                    &bsk,
                    &lwe_secret_key_after_ks,
                    polynomial_size,
                    glwe_size,
                    ciphertext_modulus,
                    delta,
                    &mut lwe_r,
                );
            r_parts[i].as_mut().copy_from_slice(lwe_r.as_ref());
        }

        let mut sum_parts: Vec<LweCiphertext<Vec<u64>>> = Vec::with_capacity(n_nibbles + 1);
        let mut ct_carry_sum = LweCiphertext::new(
            0u64,
            lwe_secret_key.lwe_dimension().to_lwe_size(),
            ciphertext_modulus,
        );
        let mut ct_tmp_sum = LweCiphertext::new(
            0u64,
            lwe_secret_key.lwe_dimension().to_lwe_size(),
            ciphertext_modulus,
        );
        encrypt_lwe_ciphertext(
            &lwe_secret_key,
            &mut ct_carry_sum,
            Plaintext(0),
            glwe_modular_std_dev,
            &mut encryption_generator,
        );
        let mut zero_ct_l = LweCiphertext::new(
            0u64,
            lwe_secret_key.lwe_dimension().to_lwe_size(),
            ciphertext_modulus,
        );
        let mut zero_ct_r = LweCiphertext::new(
            0u64,
            lwe_secret_key.lwe_dimension().to_lwe_size(),
            ciphertext_modulus,
        );
        encrypt_lwe_ciphertext(
            &lwe_secret_key,
            &mut zero_ct_l,
            Plaintext(0),
            glwe_modular_std_dev,
            &mut encryption_generator,
        );
        encrypt_lwe_ciphertext(
            &lwe_secret_key,
            &mut zero_ct_r,
            Plaintext(0),
            glwe_modular_std_dev,
            &mut encryption_generator,
        );
        for i in 0..n_nibbles + 1 {
            let x_slice = if i < n_nibbles { 
                l_parts[i].as_ref() 
            } else {
                zero_ct_l.as_ref()
            };
            
            let y_slice = if i < n_nibbles { 
                r_parts[i].as_ref() 
            } else {
                zero_ct_r.as_ref()
            };
            let carry_slice = ct_carry_sum.as_ref();
            let tmp_slice = ct_tmp_sum.as_mut();
            tmp_slice.copy_from_slice(x_slice);
            for idx in 0..tmp_slice.len() {
                tmp_slice[idx] = tmp_slice[idx].wrapping_add(y_slice[idx]);
            }
            for idx in 0..tmp_slice.len() {
                tmp_slice[idx] = tmp_slice[idx].wrapping_add(carry_slice[idx]);
            }
            pbs::pbs_carry(
                &ct_tmp_sum,
                &ksk,
                &bsk,
                &lwe_secret_key_after_ks,
                polynomial_size,
                glwe_size,
                ciphertext_modulus,
                delta,
                &mut ct_carry_sum,
            );
            let mut ct_rem = LweCiphertext::new(
                0u64,
                lwe_secret_key.lwe_dimension().to_lwe_size(),
                ciphertext_modulus,
            );            
            pbs::pbs_reminder(
                &ct_tmp_sum,
                &ksk,
                &bsk,
                &lwe_secret_key_after_ks,
                polynomial_size,
                glwe_size,
                ciphertext_modulus,
                delta,
                &mut ct_rem,
            );
            sum_parts.push(ct_rem);
        }
        let mut tmp_div2_lwe = LweCiphertext::new(
            0u64,
            lwe_secret_key.lwe_dimension().to_lwe_size(),
            ciphertext_modulus,
        );
        let mut div_carry_ct = LweCiphertext::new(
            0u64,
            lwe_secret_key.lwe_dimension().to_lwe_size(),
            ciphertext_modulus,
        );
        encrypt_lwe_ciphertext(
            &lwe_secret_key,
            &mut div_carry_ct,
            Plaintext(0),
            glwe_modular_std_dev,
            &mut encryption_generator,
        );
        for i in (0..n_nibbles + 1).rev() {
            if i == n_nibbles {
                pbs::pbs_divide_carry(
                    &sum_parts[i],
                    &ksk,
                    &bsk,
                    &lwe_secret_key_after_ks,
                    polynomial_size,
                    glwe_size,
                    ciphertext_modulus,
                    delta,
                    &mut div_carry_ct,
                );
            } else {
                pbs::pbs_divide_by_two(
                    &sum_parts[i],
                    &ksk,
                    &bsk,
                    &lwe_secret_key_after_ks,
                    polynomial_size,
                    glwe_size,
                    ciphertext_modulus,
                    delta,
                    &mut tmp_div2_lwe,
                );
                for (a, b) in tmp_div2_lwe.as_mut().iter_mut().zip(div_carry_ct.as_ref().iter()) {
                    *a = a.wrapping_add(*b);
                }
                c_parts[i].as_mut().copy_from_slice(tmp_div2_lwe.as_ref());
                pbs::pbs_divide_carry(
                    &sum_parts[i],
                    &ksk,
                    &bsk,
                    &lwe_secret_key_after_ks,
                    polynomial_size,
                    glwe_size,
                    ciphertext_modulus,
                    delta,
                    &mut div_carry_ct,
                );
            }
        }

        let mut ggsw_lists_c: Vec<GgswCiphertextList<Vec<u64>>> = Vec::with_capacity(log_n);
        for (i, ct) in c_parts.iter().enumerate()  {
            let bits = if i == n_nibbles - 1 { n_rems } else { 4 };
            let ggsw_list_out = cbs_int_6to4::cbs_int_6to4(
            ct,
            &ksk,
            &bsk,
            &auto_keys,
            &ss_key,
            bits,
            );
            ggsw_lists_c.push(ggsw_list_out);
        }
        ggsw_list_c = concat_ggsw_lists(&ggsw_lists_c);
    }
    let fourier_eq = {
        let mut fourier = FourierGgswCiphertext::new(glwe_size, polynomial_size, cbs_base_log, cbs_level);
        convert_standard_ggsw_ciphertext_to_fourier(&eq, &mut fourier);
        fourier
    };
    let error_val: u64 = 0xFu64;
    let mut ct_error = LweCiphertext::new(
        0u64,
        lwe_secret_key.lwe_dimension().to_lwe_size(),
        ciphertext_modulus,
    );
    encrypt_lwe_ciphertext(
        &lwe_secret_key,
        &mut ct_error,
        Plaintext(error_val * delta),
        glwe_modular_std_dev,
        &mut encryption_generator,
    );
    let mut sol_parts: Vec<LweCiphertext<Vec<u64>>> = Vec::with_capacity(i_parts.len());
    for i in 0..i_parts.len() {
        let mut glwe_out = GlweCiphertext::new(0u64, glwe_size, polynomial_size, ciphertext_modulus);
        let mut glwe_i = GlweCiphertext::new(0u64, glwe_size, polynomial_size, ciphertext_modulus);
        convert_lwe_to_glwe_by_trace_with_preprocessing(&ct_error, &mut glwe_out, &auto_keys);
        convert_lwe_to_glwe_by_trace_with_preprocessing(&i_parts[i], &mut glwe_i, &auto_keys);
        cmux_assign(&mut glwe_out, &mut glwe_i, &fourier_eq);
        let mut lwe_out = LweCiphertext::new(
            0u64,
            lwe_secret_key.lwe_dimension().to_lwe_size(),
            ciphertext_modulus,
        );
        extract_lwe_sample_from_glwe_ciphertext(&glwe_out, &mut lwe_out, MonomialDegree(0));
        sol_parts.push(lwe_out);
    }
    sol_parts
}