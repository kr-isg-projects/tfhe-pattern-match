use tfhe::core_crypto::fft_impl::fft64::c64;
use tfhe::core_crypto::prelude::*;
use tfhe::core_crypto::entities::FourierLweBootstrapKey;
use refined_tfhe_lhe::generate_accumulator;
use aligned_vec::ABox;

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
    keyswitch_lwe_ciphertext(
        ksk,
        lwe_in,
        &mut lwe_ks,
    );
    let mut lut = vec![0u64; 32];
    for i in 0..32 {
        lut[i] = i as u64;
    }
    let accumulator = generate_accumulator(
        polynomial_size,
        glwe_size,
        32,
        ciphertext_modulus,
        delta,
        |i| lut[i as usize],
    );
    programmable_bootstrap_lwe_ciphertext(
        &lwe_ks,
        lwe_out,
        &accumulator,
        bsk,
    );
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
    keyswitch_lwe_ciphertext(
        ksk,
        lwe_in,
        &mut lwe_ks,
    );
    let mut lut = vec![0u64; 32];
    for i in 0..32 {
        let upper = (i >> 4) & 0b1;
        lut[i as usize] = match upper {
            0b0 => 0,
            0b1 => 1,
            _ => unreachable!(),
        };
    }
    let accumulator = generate_accumulator(
        polynomial_size, 
        glwe_size,
        32, 
        ciphertext_modulus, 
        delta,
        |i| lut[i as usize], 
    );
    programmable_bootstrap_lwe_ciphertext(
        &lwe_ks,
        lwe_out_carry,
        &accumulator,
        bsk,
    );
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
    keyswitch_lwe_ciphertext( 
        ksk,
        lwe_in,
        &mut lwe_ks,
    );
    let mut lut = vec![0u64; 32];
    for i in 0..32 {
        lut[i as usize] = match i {
            0..=30 => 0,        
            31 => 63,   
            _ => unreachable!(),
        };
    }
    let accumulator = generate_accumulator(
        polynomial_size,
        glwe_size,
        32,
        ciphertext_modulus,
        delta,
        |i| lut[i as usize], 
    );
    programmable_bootstrap_lwe_ciphertext(
        &lwe_ks,
        lwe_out_carry,
        &accumulator,
        bsk,
    );
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
    keyswitch_lwe_ciphertext(
        ksk,
        lwe_in,
        &mut lwe_ks,
    );
    let mut lut = vec![0u64; 32];
    for i in 0..32 {
        lut[i] = match i {    
            0..=15 => i as u64,
            16..=31 => i as u64 - 16,
            _ => unreachable!(),
        };
    }
    let accumulator = generate_accumulator(
        polynomial_size,
        glwe_size,
        32,
        ciphertext_modulus,
        delta,
        |i| lut[i as usize], 
    );
    programmable_bootstrap_lwe_ciphertext(
        &lwe_ks,
        lwe_out_carry,
        &accumulator,
        bsk,
    );
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
    keyswitch_lwe_ciphertext(
        ksk,
        lwe_in,
        &mut lwe_ks,
    );
    let mut lut = vec![0u64; 32];
    for i in 0..32 {
        lut[i] = match i {
            0 => 0,
            1..=31 => 1,
            _ => unreachable!(),
        };
    }
    let accumulator = generate_accumulator(
        polynomial_size,
        glwe_size,
        32,
        ciphertext_modulus,
        delta,
        |i| lut[i as usize],
    );
    programmable_bootstrap_lwe_ciphertext(
        &lwe_ks,
        lwe_out,
        &accumulator,
        bsk,
    );
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
    keyswitch_lwe_ciphertext(
        ksk,
        lwe_in,
        &mut lwe_ks,
    );
    let mut lut = vec![0u64; 32];
    for i in 0..32 {
        lut[i] = match i {
            0 => 0,
            1..=4 => 1,
            5..=15 => 0,
            16..=31 => 1,
            _ => unreachable!(),
        };
    }
    let accumulator = generate_accumulator(
        polynomial_size,
        glwe_size,
        32,
        ciphertext_modulus,
        delta,
        |i| lut[i as usize],
    );
    programmable_bootstrap_lwe_ciphertext(
        &lwe_ks,
        lwe_out,
        &accumulator,
        bsk,
    );
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
    keyswitch_lwe_ciphertext(
        ksk,
        lwe_in,
        &mut lwe_ks,
    );
    let mut lut = vec![0u64; 32];
    for i in 0..32 {
        lut[i] = match i {
            0 => 0,
            2 => 1,
            4 => 1,
            28 => 1,
            29 => 0,
            30 => 63,
            31 => 0,
            _ => 0,
        };
    }
    let accumulator = generate_accumulator(
        polynomial_size,
        glwe_size,
        32,
        ciphertext_modulus,
        delta,
        |i| lut[i as usize],
    );
    programmable_bootstrap_lwe_ciphertext(
        &lwe_ks,
        lwe_out,
        &accumulator,
        bsk,
    );
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
    keyswitch_lwe_ciphertext(
        ksk,
        lwe_in,
        &mut lwe_ks,
    );
    let mut lut = vec![0u64; 32];
    for i in 0..32 {
        lut[i] = match i {
            0 => 0,  
            1 => 0,  
            2 => 1,
            3 => 1,  
            4 => 2,  
            5 => 2,  
            6 => 3,  
            7 => 3,  
            8 => 4,  
            9 => 4,  
            10 => 5,  
            11 => 5,  
            12 => 6,  
            13 => 6,  
            14 => 7,  
            15 => 7,
            _ => 0,
        };
    }
    let accumulator = generate_accumulator(
        polynomial_size,
        glwe_size,
        32,
        ciphertext_modulus,
        delta,
        |i| lut[i as usize],
    );
    programmable_bootstrap_lwe_ciphertext(
        &lwe_ks,
        lwe_out,
        &accumulator,
        bsk,
    );
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
    keyswitch_lwe_ciphertext(
        ksk,
        lwe_in,
        &mut lwe_ks,
    );
    let mut lut = vec![0u64; 32];
    for i in 0..32 {
        lut[i] = match i {
            0 => 0,  
            1 => 8,
            2 => 0,  
            3 => 8,  
            4 => 0,  
            5 => 8,  
            6 => 0,  
            7 => 8,  
            8 => 0,  
            9 => 8,  
            10 => 0,  
            11 => 8,  
            12 => 0,  
            13 => 8,  
            14 => 0,  
            15 => 8,  
            _ => 0,
        };
    }
    let accumulator = generate_accumulator(
        polynomial_size,
        glwe_size,
        32,
        ciphertext_modulus,
        delta,
        |i| lut[i as usize],
    );
    programmable_bootstrap_lwe_ciphertext(
        &lwe_ks,
        lwe_out,
        &accumulator,
        bsk,
    );
}