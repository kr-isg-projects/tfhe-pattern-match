use crate::pbs;
use aligned_vec::ABox;
use rayon::join;
use tfhe::core_crypto::algorithms::encrypt_lwe_ciphertext;
use tfhe::core_crypto::entities::FourierLweBootstrapKey;
use tfhe::core_crypto::fft_impl::fft64::c64;
use tfhe::core_crypto::prelude::*;

pub struct ArithmeticCtx<'a> {
    pub lwe_secret_key: &'a LweSecretKey<Vec<u64>>,
    pub lwe_secret_key_after_ks: &'a LweSecretKey<Vec<u64>>,
    pub ksk: &'a LweKeyswitchKey<Vec<u64>>,
    pub bsk: &'a FourierLweBootstrapKey<ABox<[c64]>>,
    pub polynomial_size: PolynomialSize,
    pub glwe_size: GlweSize,
    pub ciphertext_modulus: CiphertextModulus<u64>,
    pub delta: u64,
}

pub fn add_one_lwe_vec(
    c_parts: &[LweCiphertext<Vec<u64>>],
    noise_distribution: Gaussian<f64>,
    encryption_generator: &mut EncryptionRandomGenerator<DefaultRandomGenerator>,
    ctx: &ArithmeticCtx<'_>,
) -> Vec<LweCiphertext<Vec<u64>>> {
    let mut c_parts_p1 = c_parts.to_vec();
    let mut ct_carry_p = LweCiphertext::new(
        0u64,
        ctx.lwe_secret_key.lwe_dimension().to_lwe_size(),
        ctx.ciphertext_modulus,
    );
    encrypt_lwe_ciphertext(
        ctx.lwe_secret_key,
        &mut ct_carry_p,
        Plaintext(ctx.delta),
        noise_distribution,
        encryption_generator,
    );

    let mut ct_tmp_p = LweCiphertext::new(
        0u64,
        ctx.lwe_secret_key.lwe_dimension().to_lwe_size(),
        ctx.ciphertext_modulus,
    );
    for ct in c_parts_p1.iter_mut() {
        ct_tmp_p.as_mut().copy_from_slice(ct.as_ref());
        for (x, carry) in ct_tmp_p.as_mut().iter_mut().zip(ct_carry_p.as_ref().iter()) {
            *x = x.wrapping_add(*carry);
        }
        pbs::pbs_carry(
            &ct_tmp_p,
            ctx.ksk,
            ctx.bsk,
            ctx.lwe_secret_key_after_ks,
            ctx.polynomial_size,
            ctx.glwe_size,
            ctx.ciphertext_modulus,
            ctx.delta,
            &mut ct_carry_p,
        );
        pbs::pbs_reminder(
            &ct_tmp_p,
            ctx.ksk,
            ctx.bsk,
            ctx.lwe_secret_key_after_ks,
            ctx.polynomial_size,
            ctx.glwe_size,
            ctx.ciphertext_modulus,
            ctx.delta,
            ct,
        );
    }
    c_parts_p1
}

pub fn sub_one_lwe_vec(
    c_parts: &[LweCiphertext<Vec<u64>>],
    noise_distribution: Gaussian<f64>,
    encryption_generator: &mut EncryptionRandomGenerator<DefaultRandomGenerator>,
    ctx: &ArithmeticCtx<'_>,
) -> Vec<LweCiphertext<Vec<u64>>> {
    let mut c_parts_m1 = c_parts.to_vec();
    let mut ct_carry_m = LweCiphertext::new(
        0u64,
        ctx.lwe_secret_key.lwe_dimension().to_lwe_size(),
        ctx.ciphertext_modulus,
    );
    encrypt_lwe_ciphertext(
        ctx.lwe_secret_key,
        &mut ct_carry_m,
        Plaintext(ctx.delta),
        noise_distribution,
        encryption_generator,
    );

    let mut ct_tmp_m = LweCiphertext::new(
        0u64,
        ctx.lwe_secret_key.lwe_dimension().to_lwe_size(),
        ctx.ciphertext_modulus,
    );
    for ct in c_parts_m1.iter_mut() {
        ct_tmp_m.as_mut().copy_from_slice(ct.as_ref());
        for (x, carry) in ct_tmp_m.as_mut().iter_mut().zip(ct_carry_m.as_ref().iter()) {
            *x = x.wrapping_sub(*carry);
        }
        pbs::pbs_carry_sub(
            &ct_tmp_m,
            ctx.ksk,
            ctx.bsk,
            ctx.lwe_secret_key_after_ks,
            ctx.polynomial_size,
            ctx.glwe_size,
            ctx.ciphertext_modulus,
            ctx.delta,
            &mut ct_carry_m,
        );
        pbs::pbs_reminder(
            &ct_tmp_m,
            ctx.ksk,
            ctx.bsk,
            ctx.lwe_secret_key_after_ks,
            ctx.polynomial_size,
            ctx.glwe_size,
            ctx.ciphertext_modulus,
            ctx.delta,
            ct,
        );
    }
    c_parts_m1
}

pub fn add_lwe_vec(
    l_parts: &[LweCiphertext<Vec<u64>>],
    r_parts: &[LweCiphertext<Vec<u64>>],
    noise_distribution: Gaussian<f64>,
    encryption_generator: &mut EncryptionRandomGenerator<DefaultRandomGenerator>,
    ctx: &ArithmeticCtx<'_>,
) -> Vec<LweCiphertext<Vec<u64>>> {
    let n_parts = l_parts.len();
    let mut sum_parts: Vec<LweCiphertext<Vec<u64>>> = Vec::with_capacity(n_parts + 1);
    let mut ct_carry_sum = LweCiphertext::new(
        0u64,
        ctx.lwe_secret_key.lwe_dimension().to_lwe_size(),
        ctx.ciphertext_modulus,
    );
    let mut ct_tmp_sum = LweCiphertext::new(
        0u64,
        ctx.lwe_secret_key.lwe_dimension().to_lwe_size(),
        ctx.ciphertext_modulus,
    );
    encrypt_lwe_ciphertext(
        ctx.lwe_secret_key,
        &mut ct_carry_sum,
        Plaintext(0),
        noise_distribution,
        encryption_generator,
    );
    let mut zero_ct_l = LweCiphertext::new(
        0u64,
        ctx.lwe_secret_key.lwe_dimension().to_lwe_size(),
        ctx.ciphertext_modulus,
    );
    let mut zero_ct_r = LweCiphertext::new(
        0u64,
        ctx.lwe_secret_key.lwe_dimension().to_lwe_size(),
        ctx.ciphertext_modulus,
    );
    encrypt_lwe_ciphertext(
        ctx.lwe_secret_key,
        &mut zero_ct_l,
        Plaintext(0),
        noise_distribution,
        encryption_generator,
    );
    encrypt_lwe_ciphertext(
        ctx.lwe_secret_key,
        &mut zero_ct_r,
        Plaintext(0),
        noise_distribution,
        encryption_generator,
    );
    for i in 0..n_parts + 1 {
        let x_slice = if i < n_parts {
            l_parts[i].as_ref()
        } else {
            zero_ct_l.as_ref()
        };

        let y_slice = if i < n_parts {
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
        let mut ct_rem = LweCiphertext::new(
            0u64,
            ctx.lwe_secret_key.lwe_dimension().to_lwe_size(),
            ctx.ciphertext_modulus,
        );
        join(
            || {
                pbs::pbs_carry(
                    &ct_tmp_sum,
                    ctx.ksk,
                    ctx.bsk,
                    ctx.lwe_secret_key_after_ks,
                    ctx.polynomial_size,
                    ctx.glwe_size,
                    ctx.ciphertext_modulus,
                    ctx.delta,
                    &mut ct_carry_sum,
                );
            },
            || {
                pbs::pbs_reminder(
                    &ct_tmp_sum,
                    ctx.ksk,
                    ctx.bsk,
                    ctx.lwe_secret_key_after_ks,
                    ctx.polynomial_size,
                    ctx.glwe_size,
                    ctx.ciphertext_modulus,
                    ctx.delta,
                    &mut ct_rem,
                );
            },
        );
        sum_parts.push(ct_rem);
    }
    sum_parts
}

/// Add a known public plaintext `offset` to an encrypted radix-4 value
/// (`l_parts`, LSB-first digits), specialized for the case where the right
/// operand is cleartext rather than an encrypted vector.
///
/// Unlike `add_lwe_vec` (encrypted + encrypted, needed e.g. for `l + r` in
/// the binary search), this never encrypts `offset`: every digit below its
/// lowest non-zero radix-4 digit is left untouched with no PBS at all (the
/// offset's digit and the incoming carry are both 0 there, so the sum digit
/// equals the input digit), and from that digit onward the known plaintext
/// digit is folded in via a noiseless trivial encryption (zero mask, body =
/// digit * delta) instead of a real, randomness-consuming `encrypt_lwe_ciphertext`
/// call. The subsequent carry propagation (`pbs_carry` / `pbs_reminder`) is
/// identical to `add_lwe_vec`.
pub fn add_plain_offset_lwe_vec(
    l_parts: &[LweCiphertext<Vec<u64>>],
    offset: u64,
    ctx: &ArithmeticCtx<'_>,
) -> Vec<LweCiphertext<Vec<u64>>> {
    let n_parts = l_parts.len();
    let lwe_size = ctx.lwe_secret_key.lwe_dimension().to_lwe_size();

    // Lowest radix-4 digit index of `offset` that is non-zero. Below it,
    // both the offset digit and the incoming carry are 0 (ripple-carry
    // starts at 0), so the sum digit is just the input digit, unchanged.
    let skip = if offset == 0 {
        n_parts
    } else {
        ((offset.trailing_zeros() / 2) as usize).min(n_parts)
    };

    let mut sum_parts: Vec<LweCiphertext<Vec<u64>>> = Vec::with_capacity(n_parts + 1);
    sum_parts.extend(l_parts[..skip].iter().cloned());

    let mut ct_carry_sum =
        allocate_and_trivially_encrypt_new_lwe_ciphertext(lwe_size, Plaintext(0u64), ctx.ciphertext_modulus);

    for i in skip..=n_parts {
        let offset_digit = if i < n_parts { (offset >> (2 * i)) & 0x3 } else { 0 };

        let mut ct_tmp_sum = if i < n_parts {
            l_parts[i].clone()
        } else {
            allocate_and_trivially_encrypt_new_lwe_ciphertext(
                lwe_size,
                Plaintext(0u64),
                ctx.ciphertext_modulus,
            )
        };

        let offset_ct = allocate_and_trivially_encrypt_new_lwe_ciphertext(
            lwe_size,
            Plaintext(offset_digit * ctx.delta),
            ctx.ciphertext_modulus,
        );
        for (x, y) in ct_tmp_sum.as_mut().iter_mut().zip(offset_ct.as_ref().iter()) {
            *x = x.wrapping_add(*y);
        }
        for (x, c) in ct_tmp_sum
            .as_mut()
            .iter_mut()
            .zip(ct_carry_sum.as_ref().iter())
        {
            *x = x.wrapping_add(*c);
        }

        let mut ct_rem = LweCiphertext::new(0u64, lwe_size, ctx.ciphertext_modulus);
        let mut next_carry = LweCiphertext::new(0u64, lwe_size, ctx.ciphertext_modulus);
        join(
            || {
                pbs::pbs_carry(
                    &ct_tmp_sum,
                    ctx.ksk,
                    ctx.bsk,
                    ctx.lwe_secret_key_after_ks,
                    ctx.polynomial_size,
                    ctx.glwe_size,
                    ctx.ciphertext_modulus,
                    ctx.delta,
                    &mut next_carry,
                );
            },
            || {
                pbs::pbs_reminder(
                    &ct_tmp_sum,
                    ctx.ksk,
                    ctx.bsk,
                    ctx.lwe_secret_key_after_ks,
                    ctx.polynomial_size,
                    ctx.glwe_size,
                    ctx.ciphertext_modulus,
                    ctx.delta,
                    &mut ct_rem,
                );
            },
        );
        ct_carry_sum = next_carry;
        sum_parts.push(ct_rem);
    }

    sum_parts
}

pub fn div2_lwe_vec(
    sum_parts: &[LweCiphertext<Vec<u64>>],
    n_parts: usize,
    noise_distribution: Gaussian<f64>,
    encryption_generator: &mut EncryptionRandomGenerator<DefaultRandomGenerator>,
    ctx: &ArithmeticCtx<'_>,
) -> Vec<LweCiphertext<Vec<u64>>> {
    let mut out_parts: Vec<LweCiphertext<Vec<u64>>> = sum_parts[..n_parts].to_vec();
    let mut tmp_div2_lwe = LweCiphertext::new(
        0u64,
        ctx.lwe_secret_key.lwe_dimension().to_lwe_size(),
        ctx.ciphertext_modulus,
    );
    let mut div_carry_ct = LweCiphertext::new(
        0u64,
        ctx.lwe_secret_key.lwe_dimension().to_lwe_size(),
        ctx.ciphertext_modulus,
    );
    encrypt_lwe_ciphertext(
        ctx.lwe_secret_key,
        &mut div_carry_ct,
        Plaintext(0),
        noise_distribution,
        encryption_generator,
    );
    for i in (0..n_parts + 1).rev() {
        if i == n_parts {
            pbs::pbs_divide_carry(
                &sum_parts[i],
                ctx.ksk,
                ctx.bsk,
                ctx.lwe_secret_key_after_ks,
                ctx.polynomial_size,
                ctx.glwe_size,
                ctx.ciphertext_modulus,
                ctx.delta,
                &mut div_carry_ct,
            );
        } else {
            let incoming_carry = div_carry_ct.clone();
            let mut next_div_carry = LweCiphertext::new(
                0u64,
                ctx.lwe_secret_key.lwe_dimension().to_lwe_size(),
                ctx.ciphertext_modulus,
            );
            join(
                || {
                    pbs::pbs_divide_by_two(
                        &sum_parts[i],
                        ctx.ksk,
                        ctx.bsk,
                        ctx.lwe_secret_key_after_ks,
                        ctx.polynomial_size,
                        ctx.glwe_size,
                        ctx.ciphertext_modulus,
                        ctx.delta,
                        &mut tmp_div2_lwe,
                    );
                },
                || {
                    pbs::pbs_divide_carry(
                        &sum_parts[i],
                        ctx.ksk,
                        ctx.bsk,
                        ctx.lwe_secret_key_after_ks,
                        ctx.polynomial_size,
                        ctx.glwe_size,
                        ctx.ciphertext_modulus,
                        ctx.delta,
                        &mut next_div_carry,
                    );
                },
            );
            for (a, b) in tmp_div2_lwe
                .as_mut()
                .iter_mut()
                .zip(incoming_carry.as_ref().iter())
            {
                *a = a.wrapping_add(*b);
            }
            out_parts[i].as_mut().copy_from_slice(tmp_div2_lwe.as_ref());
            div_carry_ct = next_div_carry;
        }
    }
    out_parts
}
