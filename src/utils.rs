use tfhe::core_crypto::algorithms::encrypt_lwe_ciphertext;
use tfhe::core_crypto::prelude::*;

pub fn concat_ggsw_lists(lists: &[GgswCiphertextList<Vec<u64>>]) -> GgswCiphertextList<Vec<u64>> {
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

pub fn encrypt_small_int(
    value: u64,
    lwe_secret_key: &LweSecretKey<Vec<u64>>,
    ciphertext_modulus: CiphertextModulus<u64>,
    noise_distribution: Gaussian<f64>,
    encryption_generator: &mut EncryptionRandomGenerator<DefaultRandomGenerator>,
    delta: u64,
) -> LweCiphertext<Vec<u64>> {
    let mut ct = LweCiphertext::new(
        0u64,
        lwe_secret_key.lwe_dimension().to_lwe_size(),
        ciphertext_modulus,
    );
    encrypt_lwe_ciphertext(
        lwe_secret_key,
        &mut ct,
        Plaintext(value * delta),
        noise_distribution,
        encryption_generator,
    );
    ct
}

pub fn init_constant_parts_lwe(
    value: u64,
    n_parts: usize,
    lwe_secret_key: &LweSecretKey<Vec<u64>>,
    ciphertext_modulus: CiphertextModulus<u64>,
    noise_distribution: Gaussian<f64>,
    encryption_generator: &mut EncryptionRandomGenerator<DefaultRandomGenerator>,
    delta: u64,
) -> Vec<LweCiphertext<Vec<u64>>> {
    (0..n_parts)
        .map(|_| {
            encrypt_small_int(
                value,
                lwe_secret_key,
                ciphertext_modulus,
                noise_distribution,
                encryption_generator,
                delta,
            )
        })
        .collect()
}

pub fn init_radix_parts_lwe(
    value: u64,
    n_parts: usize,
    lwe_secret_key: &LweSecretKey<Vec<u64>>,
    ciphertext_modulus: CiphertextModulus<u64>,
    noise_distribution: Gaussian<f64>,
    encryption_generator: &mut EncryptionRandomGenerator<DefaultRandomGenerator>,
    delta: u64,
) -> Vec<LweCiphertext<Vec<u64>>> {
    (0..n_parts)
        .map(|i| {
            let part_val = (value >> (i * 2)) & 0x3;
            encrypt_small_int(
                part_val,
                lwe_secret_key,
                ciphertext_modulus,
                noise_distribution,
                encryption_generator,
                delta,
            )
        })
        .collect()
}
