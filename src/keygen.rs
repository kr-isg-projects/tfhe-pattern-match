use tfhe::core_crypto::prelude::*;

pub fn keygen_pbs<G: ByteRandomGenerator>(
    lwe_dimension: LweDimension,
    glwe_dimension: GlweDimension,
    polynomial_size: PolynomialSize,
    lwe_modular_std_dev: impl DispersionParameter,
    glwe_modular_std_dev: impl DispersionParameter,
    pbs_base_log: DecompositionBaseLog,
    pbs_level: DecompositionLevelCount,
    ks_base_log: DecompositionBaseLog,
    ks_level: DecompositionLevelCount,
    secret_generator: &mut SecretRandomGenerator<G>,
    encryption_generator: &mut EncryptionRandomGenerator<G>,
) -> (
    LweSecretKey<Vec<u64>>,
    GlweSecretKey<Vec<u64>>,
    LweSecretKey<Vec<u64>>,
    FourierLweBootstrapKeyOwned,
    LweKeyswitchKey<Vec<u64>>,
) {
    let small_lwe_secret_key = LweSecretKey::generate_new_binary(lwe_dimension, secret_generator);
    let glwe_secret_key =
        GlweSecretKey::generate_new_binary(glwe_dimension, polynomial_size, secret_generator);
    let lwe_secret_key = glwe_secret_key.clone().into_lwe_secret_key();
    let lwe_secret_key_after_ks = small_lwe_secret_key;

    let bootstrap_key = allocate_and_generate_new_lwe_bootstrap_key(
        &lwe_secret_key_after_ks,
        &glwe_secret_key,
        pbs_base_log,
        pbs_level,
        Gaussian::from_dispersion_parameter(glwe_modular_std_dev, 0.0),
        CiphertextModulus::<u64>::new_native(),
        encryption_generator,
    );

    let mut fourier_bsk = FourierLweBootstrapKey::new(
        bootstrap_key.input_lwe_dimension(),
        bootstrap_key.glwe_size(),
        bootstrap_key.polynomial_size(),
        bootstrap_key.decomposition_base_log(),
        bootstrap_key.decomposition_level_count(),
    );
    convert_standard_lwe_bootstrap_key_to_fourier(&bootstrap_key, &mut fourier_bsk);

    let ksk = allocate_and_generate_new_lwe_keyswitch_key(
        &lwe_secret_key,
        &lwe_secret_key_after_ks,
        ks_base_log,
        ks_level,
        Gaussian::from_dispersion_parameter(lwe_modular_std_dev, 0.0),
        CiphertextModulus::<u64>::new_native(),
        encryption_generator,
    );

    (
        lwe_secret_key,
        glwe_secret_key,
        lwe_secret_key_after_ks,
        fourier_bsk,
        ksk,
    )
}
