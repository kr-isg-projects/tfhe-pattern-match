use aligned_vec::ABox;
use tfhe::core_crypto::entities::FourierLweBootstrapKey;
use tfhe::core_crypto::fft_impl::fft64::c64;
use tfhe::core_crypto::prelude::*;

const SIGNED_LUT_SIZE: usize = 8;

pub fn generate_accumulator<F>(
    polynomial_size: PolynomialSize,
    glwe_size: GlweSize,
    message_modulus: usize,
    ciphertext_modulus: CiphertextModulus<u64>,
    delta: u64,
    f: F,
) -> GlweCiphertextOwned<u64>
where
    F: Fn(u64) -> u64,
{
    let box_size = polynomial_size.0 / message_modulus;
    let mut accumulator = vec![0u64; polynomial_size.0];

    for i in 0..message_modulus {
        let index = i * box_size;
        accumulator[index..index + box_size]
            .iter_mut()
            .for_each(|value| *value = f(i as u64) * delta);
    }

    let half_box_size = box_size / 2;
    if ciphertext_modulus.is_compatible_with_native_modulus() {
        for value in &mut accumulator[0..half_box_size] {
            *value = value.wrapping_neg();
        }
    } else {
        let modulus = ciphertext_modulus.get_custom_modulus() as u64;
        for value in &mut accumulator[0..half_box_size] {
            *value = value.wrapping_neg_custom_mod(modulus);
        }
    }

    accumulator.rotate_left(half_box_size);
    let accumulator_plaintext = PlaintextList::from_container(accumulator);
    allocate_and_trivially_encrypt_new_glwe_ciphertext(
        glwe_size,
        &accumulator_plaintext,
        ciphertext_modulus,
    )
}

pub fn pbs_identity(
    lwe_in: &LweCiphertext<Vec<u64>>,
    ksk: &LweKeyswitchKey<Vec<u64>>,
    bsk: &FourierLweBootstrapKey<ABox<[c64]>>,
    lwe_secret_key_after_ks: &LweSecretKey<Vec<u64>>,
    polynomial_size: PolynomialSize,
    glwe_size: GlweSize,
    ciphertext_modulus: CiphertextModulus<u64>,
    delta: u64,
    lwe_out: &mut LweCiphertext<Vec<u64>>,
) {
    let mut lwe_ks = LweCiphertext::new(
        0u64,
        lwe_secret_key_after_ks.lwe_dimension().to_lwe_size(),
        ciphertext_modulus,
    );
    keyswitch_lwe_ciphertext(ksk, lwe_in, &mut lwe_ks);
    let mut lut = vec![0u64; SIGNED_LUT_SIZE];
    for i in 0..SIGNED_LUT_SIZE {
        lut[i] = i as u64;
    }
    let accumulator = generate_accumulator(
        polynomial_size,
        glwe_size,
        SIGNED_LUT_SIZE,
        ciphertext_modulus,
        delta,
        |i| lut[i as usize],
    );
    programmable_bootstrap_lwe_ciphertext(&lwe_ks, lwe_out, &accumulator, bsk);
}

pub fn pbs_carry(
    lwe_in: &LweCiphertext<Vec<u64>>,
    ksk: &LweKeyswitchKey<Vec<u64>>,
    bsk: &FourierLweBootstrapKey<ABox<[c64]>>,
    lwe_secret_key_after_ks: &LweSecretKey<Vec<u64>>,
    polynomial_size: PolynomialSize,
    glwe_size: GlweSize,
    ciphertext_modulus: CiphertextModulus<u64>,
    delta: u64,
    lwe_out_carry: &mut LweCiphertext<Vec<u64>>,
) {
    let mut lwe_ks = LweCiphertext::new(
        0u64,
        lwe_secret_key_after_ks.lwe_dimension().to_lwe_size(),
        ciphertext_modulus,
    );
    keyswitch_lwe_ciphertext(ksk, lwe_in, &mut lwe_ks);
    let mut lut = vec![0u64; SIGNED_LUT_SIZE];
    for i in 0..SIGNED_LUT_SIZE {
        let upper = (i >> 2) & 0b1;
        lut[i as usize] = match upper {
            0b0 => 0,
            0b1 => 1,
            _ => unreachable!(),
        };
    }
    let accumulator = generate_accumulator(
        polynomial_size,
        glwe_size,
        SIGNED_LUT_SIZE,
        ciphertext_modulus,
        delta,
        |i| lut[i as usize],
    );
    programmable_bootstrap_lwe_ciphertext(&lwe_ks, lwe_out_carry, &accumulator, bsk);
}

pub fn pbs_carry_sub(
    lwe_in: &LweCiphertext<Vec<u64>>,
    ksk: &LweKeyswitchKey<Vec<u64>>,
    bsk: &FourierLweBootstrapKey<ABox<[c64]>>,
    lwe_secret_key_after_ks: &LweSecretKey<Vec<u64>>,
    polynomial_size: PolynomialSize,
    glwe_size: GlweSize,
    ciphertext_modulus: CiphertextModulus<u64>,
    delta: u64,
    lwe_out_carry: &mut LweCiphertext<Vec<u64>>,
) {
    let mut lwe_ks = LweCiphertext::new(
        0u64,
        lwe_secret_key_after_ks.lwe_dimension().to_lwe_size(),
        ciphertext_modulus,
    );
    keyswitch_lwe_ciphertext(ksk, lwe_in, &mut lwe_ks);
    let mut lut = vec![0u64; SIGNED_LUT_SIZE];
    for i in 0..SIGNED_LUT_SIZE {
        lut[i as usize] = match i {
            0..=6 => 0,
            7 => 3,
            _ => unreachable!(),
        };
    }
    let accumulator = generate_accumulator(
        polynomial_size,
        glwe_size,
        SIGNED_LUT_SIZE,
        ciphertext_modulus,
        delta,
        |i| lut[i as usize],
    );
    programmable_bootstrap_lwe_ciphertext(&lwe_ks, lwe_out_carry, &accumulator, bsk);
}

pub fn pbs_reminder(
    lwe_in: &LweCiphertext<Vec<u64>>,
    ksk: &LweKeyswitchKey<Vec<u64>>,
    bsk: &FourierLweBootstrapKey<ABox<[c64]>>,
    lwe_secret_key_after_ks: &LweSecretKey<Vec<u64>>,
    polynomial_size: PolynomialSize,
    glwe_size: GlweSize,
    ciphertext_modulus: CiphertextModulus<u64>,
    delta: u64,
    lwe_out_carry: &mut LweCiphertext<Vec<u64>>,
) {
    let mut lwe_ks = LweCiphertext::new(
        0u64,
        lwe_secret_key_after_ks.lwe_dimension().to_lwe_size(),
        ciphertext_modulus,
    );
    keyswitch_lwe_ciphertext(ksk, lwe_in, &mut lwe_ks);
    let mut lut = vec![0u64; SIGNED_LUT_SIZE];
    for i in 0..SIGNED_LUT_SIZE {
        lut[i] = match i {
            0..=3 => i as u64,
            4..=7 => i as u64 - 4,
            _ => unreachable!(),
        };
    }
    let accumulator = generate_accumulator(
        polynomial_size,
        glwe_size,
        SIGNED_LUT_SIZE,
        ciphertext_modulus,
        delta,
        |i| lut[i as usize],
    );
    programmable_bootstrap_lwe_ciphertext(&lwe_ks, lwe_out_carry, &accumulator, bsk);
}

pub fn pbs_sgn(
    lwe_in: &LweCiphertext<Vec<u64>>,
    ksk: &LweKeyswitchKey<Vec<u64>>,
    bsk: &FourierLweBootstrapKey<ABox<[c64]>>,
    lwe_secret_key_after_ks: &LweSecretKey<Vec<u64>>,
    polynomial_size: PolynomialSize,
    glwe_size: GlweSize,
    ciphertext_modulus: CiphertextModulus<u64>,
    delta: u64,
    lwe_out: &mut LweCiphertext<Vec<u64>>,
) {
    let mut lwe_ks = LweCiphertext::new(
        0u64,
        lwe_secret_key_after_ks.lwe_dimension().to_lwe_size(),
        ciphertext_modulus,
    );
    keyswitch_lwe_ciphertext(ksk, lwe_in, &mut lwe_ks);
    let mut lut = vec![0u64; SIGNED_LUT_SIZE];
    for i in 0..SIGNED_LUT_SIZE {
        lut[i] = match i {
            0 => 0,
            1..=7 => 1,
            _ => unreachable!(),
        };
    }
    let accumulator = generate_accumulator(
        polynomial_size,
        glwe_size,
        SIGNED_LUT_SIZE,
        ciphertext_modulus,
        delta,
        |i| lut[i as usize],
    );
    programmable_bootstrap_lwe_ciphertext(&lwe_ks, lwe_out, &accumulator, bsk);
}

pub fn pbs_cmp_pm1_zero(
    lwe_in: &LweCiphertext<Vec<u64>>,
    ksk: &LweKeyswitchKey<Vec<u64>>,
    bsk: &FourierLweBootstrapKey<ABox<[c64]>>,
    lwe_secret_key_after_ks: &LweSecretKey<Vec<u64>>,
    polynomial_size: PolynomialSize,
    glwe_size: GlweSize,
    ciphertext_modulus: CiphertextModulus<u64>,
    delta: u64,
    lwe_out: &mut LweCiphertext<Vec<u64>>,
) {
    let mut lwe_ks = LweCiphertext::new(
        0u64,
        lwe_secret_key_after_ks.lwe_dimension().to_lwe_size(),
        ciphertext_modulus,
    );
    keyswitch_lwe_ciphertext(ksk, lwe_in, &mut lwe_ks);
    let mut lut = vec![0u64; SIGNED_LUT_SIZE];
    for i in 0..SIGNED_LUT_SIZE {
        lut[i] = match i {
            0 => 0,
            1..=7 => 1,
            _ => unreachable!(),
        };
    }
    let accumulator = generate_accumulator(
        polynomial_size,
        glwe_size,
        SIGNED_LUT_SIZE,
        ciphertext_modulus,
        delta,
        |i| lut[i as usize],
    );
    programmable_bootstrap_lwe_ciphertext(&lwe_ks, lwe_out, &accumulator, bsk);
}

pub fn pbs_is_zero(
    lwe_in: &LweCiphertext<Vec<u64>>,
    ksk: &LweKeyswitchKey<Vec<u64>>,
    bsk: &FourierLweBootstrapKey<ABox<[c64]>>,
    lwe_secret_key_after_ks: &LweSecretKey<Vec<u64>>,
    polynomial_size: PolynomialSize,
    glwe_size: GlweSize,
    ciphertext_modulus: CiphertextModulus<u64>,
    delta: u64,
    lwe_out: &mut LweCiphertext<Vec<u64>>,
) {
    let mut lwe_ks = LweCiphertext::new(
        0u64,
        lwe_secret_key_after_ks.lwe_dimension().to_lwe_size(),
        ciphertext_modulus,
    );
    keyswitch_lwe_ciphertext(ksk, lwe_in, &mut lwe_ks);
    let mut lut = vec![0u64; SIGNED_LUT_SIZE];
    for (i, val) in lut.iter_mut().enumerate() {
        *val = if i == 0 { 1 } else { 0 };
    }
    let accumulator = generate_accumulator(
        polynomial_size,
        glwe_size,
        SIGNED_LUT_SIZE,
        ciphertext_modulus,
        delta,
        |i| lut[i as usize],
    );
    programmable_bootstrap_lwe_ciphertext(&lwe_ks, lwe_out, &accumulator, bsk);
}

pub fn pbs_or_match_wild(
    lwe_in: &LweCiphertext<Vec<u64>>,
    ksk: &LweKeyswitchKey<Vec<u64>>,
    bsk: &FourierLweBootstrapKey<ABox<[c64]>>,
    lwe_secret_key_after_ks: &LweSecretKey<Vec<u64>>,
    polynomial_size: PolynomialSize,
    glwe_size: GlweSize,
    ciphertext_modulus: CiphertextModulus<u64>,
    delta: u64,
    lwe_out: &mut LweCiphertext<Vec<u64>>,
) {
    let mut lwe_ks = LweCiphertext::new(
        0u64,
        lwe_secret_key_after_ks.lwe_dimension().to_lwe_size(),
        ciphertext_modulus,
    );
    keyswitch_lwe_ciphertext(ksk, lwe_in, &mut lwe_ks);
    let mut lut = vec![0u64; SIGNED_LUT_SIZE];
    for (i, val) in lut.iter_mut().enumerate() {
        *val = if i == 0 { 0 } else { 1 };
    }
    let accumulator = generate_accumulator(
        polynomial_size,
        glwe_size,
        SIGNED_LUT_SIZE,
        ciphertext_modulus,
        delta,
        |i| lut[i as usize],
    );
    programmable_bootstrap_lwe_ciphertext(&lwe_ks, lwe_out, &accumulator, bsk);
}

pub fn pbs_not_bit(
    lwe_in: &LweCiphertext<Vec<u64>>,
    ksk: &LweKeyswitchKey<Vec<u64>>,
    bsk: &FourierLweBootstrapKey<ABox<[c64]>>,
    lwe_secret_key_after_ks: &LweSecretKey<Vec<u64>>,
    polynomial_size: PolynomialSize,
    glwe_size: GlweSize,
    ciphertext_modulus: CiphertextModulus<u64>,
    delta: u64,
    lwe_out: &mut LweCiphertext<Vec<u64>>,
) {
    let mut lwe_ks = LweCiphertext::new(
        0u64,
        lwe_secret_key_after_ks.lwe_dimension().to_lwe_size(),
        ciphertext_modulus,
    );
    keyswitch_lwe_ciphertext(ksk, lwe_in, &mut lwe_ks);
    let mut lut = vec![0u64; SIGNED_LUT_SIZE];
    for (i, val) in lut.iter_mut().enumerate() {
        *val = match i {
            0 => 1,
            1 => 0,
            _ => 0,
        };
    }
    let accumulator = generate_accumulator(
        polynomial_size,
        glwe_size,
        SIGNED_LUT_SIZE,
        ciphertext_modulus,
        delta,
        |i| lut[i as usize],
    );
    programmable_bootstrap_lwe_ciphertext(&lwe_ks, lwe_out, &accumulator, bsk);
}

pub fn pbs_selector_pm1_from_last_char(
    lwe_in: &LweCiphertext<Vec<u64>>,
    ksk: &LweKeyswitchKey<Vec<u64>>,
    bsk: &FourierLweBootstrapKey<ABox<[c64]>>,
    lwe_secret_key_after_ks: &LweSecretKey<Vec<u64>>,
    polynomial_size: PolynomialSize,
    glwe_size: GlweSize,
    ciphertext_modulus: CiphertextModulus<u64>,
    delta: u64,
    lwe_out: &mut LweCiphertext<Vec<u64>>,
) {
    let mut lwe_ks = LweCiphertext::new(
        0u64,
        lwe_secret_key_after_ks.lwe_dimension().to_lwe_size(),
        ciphertext_modulus,
    );
    keyswitch_lwe_ciphertext(ksk, lwe_in, &mut lwe_ks);
    let mut lut = vec![0u64; SIGNED_LUT_SIZE];
    for val in lut.iter_mut() {
        // With negacyclic PBS, constant +1 on [0..7] yields +1 for non-negative
        // and -1 for negative inputs (encoded as 15 in 4-bit space).
        *val = 1;
    }
    let accumulator = generate_accumulator(
        polynomial_size,
        glwe_size,
        SIGNED_LUT_SIZE,
        ciphertext_modulus,
        delta,
        |i| lut[i as usize],
    );
    programmable_bootstrap_lwe_ciphertext(&lwe_ks, lwe_out, &accumulator, bsk);
}

pub fn pbs_selector_bit3(
    lwe_in: &LweCiphertext<Vec<u64>>,
    ksk: &LweKeyswitchKey<Vec<u64>>,
    bsk: &FourierLweBootstrapKey<ABox<[c64]>>,
    lwe_secret_key_after_ks: &LweSecretKey<Vec<u64>>,
    polynomial_size: PolynomialSize,
    glwe_size: GlweSize,
    ciphertext_modulus: CiphertextModulus<u64>,
    delta: u64,
    lwe_out: &mut LweCiphertext<Vec<u64>>,
) {
    let mut lwe_ks = LweCiphertext::new(
        0u64,
        lwe_secret_key_after_ks.lwe_dimension().to_lwe_size(),
        ciphertext_modulus,
    );
    keyswitch_lwe_ciphertext(ksk, lwe_in, &mut lwe_ks);
    let mut lut = vec![0u64; SIGNED_LUT_SIZE];
    for (i, val) in lut.iter_mut().enumerate() {
        *val = ((i >> 3) & 1) as u64;
    }
    let accumulator = generate_accumulator(
        polynomial_size,
        glwe_size,
        SIGNED_LUT_SIZE,
        ciphertext_modulus,
        delta,
        |i| lut[i as usize],
    );
    programmable_bootstrap_lwe_ciphertext(&lwe_ks, lwe_out, &accumulator, bsk);
}

pub fn pbs_sgn_with_special_words(
    lwe_in: &LweCiphertext<Vec<u64>>,
    ksk: &LweKeyswitchKey<Vec<u64>>,
    bsk: &FourierLweBootstrapKey<ABox<[c64]>>,
    lwe_secret_key_after_ks: &LweSecretKey<Vec<u64>>,
    polynomial_size: PolynomialSize,
    glwe_size: GlweSize,
    ciphertext_modulus: CiphertextModulus<u64>,
    delta: u64,
    lwe_out: &mut LweCiphertext<Vec<u64>>,
) {
    let mut lwe_ks = LweCiphertext::new(
        0u64,
        lwe_secret_key_after_ks.lwe_dimension().to_lwe_size(),
        ciphertext_modulus,
    );
    keyswitch_lwe_ciphertext(ksk, lwe_in, &mut lwe_ks);
    let mut lut = vec![0u64; SIGNED_LUT_SIZE];
    for i in 0..SIGNED_LUT_SIZE {
        lut[i] = match i {
            0 => 0,
            1..=4 => 1,
            5..=7 => 0,
            _ => unreachable!(),
        };
    }
    let accumulator = generate_accumulator(
        polynomial_size,
        glwe_size,
        SIGNED_LUT_SIZE,
        ciphertext_modulus,
        delta,
        |i| lut[i as usize],
    );
    programmable_bootstrap_lwe_ciphertext(&lwe_ks, lwe_out, &accumulator, bsk);
}

pub fn pbs_true_sgn(
    lwe_in: &LweCiphertext<Vec<u64>>,
    ksk: &LweKeyswitchKey<Vec<u64>>,
    bsk: &FourierLweBootstrapKey<ABox<[c64]>>,
    lwe_secret_key_after_ks: &LweSecretKey<Vec<u64>>,
    polynomial_size: PolynomialSize,
    glwe_size: GlweSize,
    ciphertext_modulus: CiphertextModulus<u64>,
    delta: u64,
    lwe_out: &mut LweCiphertext<Vec<u64>>,
) {
    let mut lwe_ks = LweCiphertext::new(
        0u64,
        lwe_secret_key_after_ks.lwe_dimension().to_lwe_size(),
        ciphertext_modulus,
    );
    keyswitch_lwe_ciphertext(ksk, lwe_in, &mut lwe_ks);
    let mut lut = vec![0u64; SIGNED_LUT_SIZE];
    for i in 0..SIGNED_LUT_SIZE {
        lut[i] = match i {
            0 => 0,
            1..=7 => 1,
            _ => 0,
        };
    }
    let accumulator = generate_accumulator(
        polynomial_size,
        glwe_size,
        SIGNED_LUT_SIZE,
        ciphertext_modulus,
        delta,
        |i| lut[i as usize],
    );
    programmable_bootstrap_lwe_ciphertext(&lwe_ks, lwe_out, &accumulator, bsk);
}

pub fn pbs_divide_by_two(
    lwe_in: &LweCiphertext<Vec<u64>>,
    ksk: &LweKeyswitchKey<Vec<u64>>,
    bsk: &FourierLweBootstrapKey<ABox<[c64]>>,
    lwe_secret_key_after_ks: &LweSecretKey<Vec<u64>>,
    polynomial_size: PolynomialSize,
    glwe_size: GlweSize,
    ciphertext_modulus: CiphertextModulus<u64>,
    delta: u64,
    lwe_out: &mut LweCiphertext<Vec<u64>>,
) {
    let mut lwe_ks = LweCiphertext::new(
        0u64,
        lwe_secret_key_after_ks.lwe_dimension().to_lwe_size(),
        ciphertext_modulus,
    );
    keyswitch_lwe_ciphertext(ksk, lwe_in, &mut lwe_ks);
    let mut lut = vec![0u64; SIGNED_LUT_SIZE];
    for i in 0..SIGNED_LUT_SIZE {
        lut[i] = match i {
            0 => 0,
            1 => 0,
            2 => 1,
            3 => 1,
            4 => 2,
            5 => 2,
            6 => 3,
            7 => 3,
            _ => 0,
        };
    }
    let accumulator = generate_accumulator(
        polynomial_size,
        glwe_size,
        SIGNED_LUT_SIZE,
        ciphertext_modulus,
        delta,
        |i| lut[i as usize],
    );
    programmable_bootstrap_lwe_ciphertext(&lwe_ks, lwe_out, &accumulator, bsk);
}

pub fn pbs_divide_carry(
    lwe_in: &LweCiphertext<Vec<u64>>,
    ksk: &LweKeyswitchKey<Vec<u64>>,
    bsk: &FourierLweBootstrapKey<ABox<[c64]>>,
    lwe_secret_key_after_ks: &LweSecretKey<Vec<u64>>,
    polynomial_size: PolynomialSize,
    glwe_size: GlweSize,
    ciphertext_modulus: CiphertextModulus<u64>,
    delta: u64,
    lwe_out: &mut LweCiphertext<Vec<u64>>,
) {
    let mut lwe_ks = LweCiphertext::new(
        0u64,
        lwe_secret_key_after_ks.lwe_dimension().to_lwe_size(),
        ciphertext_modulus,
    );
    keyswitch_lwe_ciphertext(ksk, lwe_in, &mut lwe_ks);
    let mut lut = vec![0u64; SIGNED_LUT_SIZE];
    for i in 0..SIGNED_LUT_SIZE {
        lut[i] = match i {
            0 => 0,
            1 => 2,
            2 => 0,
            3 => 2,
            4 => 0,
            5 => 2,
            6 => 0,
            7 => 2,
            _ => 0,
        };
    }
    let accumulator = generate_accumulator(
        polynomial_size,
        glwe_size,
        SIGNED_LUT_SIZE,
        ciphertext_modulus,
        delta,
        |i| lut[i as usize],
    );
    programmable_bootstrap_lwe_ciphertext(&lwe_ks, lwe_out, &accumulator, bsk);
}

pub fn pbs_mask_wild_cmp(
    lwe_in: &LweCiphertext<Vec<u64>>,
    ksk: &LweKeyswitchKey<Vec<u64>>,
    bsk: &FourierLweBootstrapKey<ABox<[c64]>>,
    lwe_secret_key_after_ks: &LweSecretKey<Vec<u64>>,
    polynomial_size: PolynomialSize,
    glwe_size: GlweSize,
    ciphertext_modulus: CiphertextModulus<u64>,
    delta: u64,
    lwe_out: &mut LweCiphertext<Vec<u64>>,
) {
    let mut lwe_ks = LweCiphertext::new(
        0u64,
        lwe_secret_key_after_ks.lwe_dimension().to_lwe_size(),
        ciphertext_modulus,
    );

    keyswitch_lwe_ciphertext(ksk, lwe_in, &mut lwe_ks);

    // Required:
    //   0  -> 0
    //   1  -> 1
    //   2  -> 0
    //   3  -> 0
    //   4  -> 0
    //   15 -> 15
    //
    // Negacyclic relation:
    //   f(x + 8) = -f(x)
    //
    // Therefore f(15)=15 is obtained from f(7)=1.
    let lut = [
        0u64, // 0
        1u64, // 1
        0u64, // 2
        0u64, // 3
        0u64, // 4
        0u64, // 5
        0u64, // 6
        1u64, // 7 -> gives f(15) = -1 = 15 mod 16
    ];

    let accumulator = generate_accumulator(
        polynomial_size,
        glwe_size,
        SIGNED_LUT_SIZE,
        ciphertext_modulus,
        delta,
        |i| lut[i as usize],
    );

    programmable_bootstrap_lwe_ciphertext(&lwe_ks, lwe_out, &accumulator, bsk);
}

pub fn pbs_eq_one(
    lwe_in: &LweCiphertext<Vec<u64>>,
    ksk: &LweKeyswitchKey<Vec<u64>>,
    bsk: &FourierLweBootstrapKey<ABox<[c64]>>,
    lwe_secret_key_after_ks: &LweSecretKey<Vec<u64>>,
    polynomial_size: PolynomialSize,
    glwe_size: GlweSize,
    ciphertext_modulus: CiphertextModulus<u64>,
    delta: u64,
    lwe_out: &mut LweCiphertext<Vec<u64>>,
) {
    let mut lwe_ks = LweCiphertext::new(
        0u64,
        lwe_secret_key_after_ks.lwe_dimension().to_lwe_size(),
        ciphertext_modulus,
    );
    keyswitch_lwe_ciphertext(ksk, lwe_in, &mut lwe_ks);
    let mut lut = vec![0u64; SIGNED_LUT_SIZE];
    for i in 0..SIGNED_LUT_SIZE {
        lut[i] = match i {
            1 => 1,
            _ => 0,
        };
    }
    let accumulator = generate_accumulator(
        polynomial_size,
        glwe_size,
        SIGNED_LUT_SIZE,
        ciphertext_modulus,
        delta,
        |i| lut[i as usize],
    );
    programmable_bootstrap_lwe_ciphertext(&lwe_ks, lwe_out, &accumulator, bsk);
}

pub fn pbs_first_nonzero(
    lwe_in: &LweCiphertext<Vec<u64>>,
    ksk: &LweKeyswitchKey<Vec<u64>>,
    bsk: &FourierLweBootstrapKey<ABox<[c64]>>,
    lwe_secret_key_after_ks: &LweSecretKey<Vec<u64>>,
    polynomial_size: PolynomialSize,
    glwe_size: GlweSize,
    ciphertext_modulus: CiphertextModulus<u64>,
    delta: u64,
    lwe_out: &mut LweCiphertext<Vec<u64>>,
) {
    let mut lwe_ks = LweCiphertext::new(
        0u64,
        lwe_secret_key_after_ks.lwe_dimension().to_lwe_size(),
        ciphertext_modulus,
    );

    keyswitch_lwe_ciphertext(ksk, lwe_in, &mut lwe_ks);

    // Input encoding:
    //
    // x = left + 3 * right
    //
    // Required direct mappings:
    //
    //   0  -> 0
    //   1  -> 1
    //   2  -> 15
    //   3  -> 1
    //   4  -> 1
    //   12 -> 15
    //   13 -> 15
    //   14 -> 1
    //   15 -> 15
    //
    // Using f(x + 8) = -f(x), the 8-entry LUT is:
    //
    //   f(0) = 0
    //   f(1) = 1
    //   f(2) = 15
    //   f(3) = 1
    //   f(4) = 1       -> f(12) = 15
    //   f(5) = 1       -> f(13) = 15
    //   f(6) = 15      -> f(14) = 1
    //   f(7) = 1       -> f(15) = 15
    let lut = [0u64, 1u64, 15u64, 1u64, 1u64, 1u64, 15u64, 1u64];

    let accumulator = generate_accumulator(
        polynomial_size,
        glwe_size,
        SIGNED_LUT_SIZE,
        ciphertext_modulus,
        delta,
        |i| lut[i as usize],
    );

    programmable_bootstrap_lwe_ciphertext(&lwe_ks, lwe_out, &accumulator, bsk);
}
