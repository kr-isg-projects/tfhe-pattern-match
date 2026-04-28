use tfhe::core_crypto::prelude::*;
use tfhe::core_crypto::fft_impl::fft64::c64;
use refined_tfhe_lhe::AutomorphKey;
use refined_tfhe_lhe::FftType;
use refined_tfhe_lhe::int_lhe_params::IntLheParam;
use aligned_vec::ABox;
use std::collections::HashMap;

pub const STD_DEV_808: StandardDev = StandardDev(2.1124945159091033e-06);
pub const STD_DEV_4096: StandardDev = StandardDev(2.168404344971009e-19);

lazy_static::lazy_static! {
    pub static ref HYBRID_BASE_64: IntLheParam<u64> = IntLheParam::new(
        LweDimension(808),
        STD_DEV_808,
        PolynomialSize(4096),
        GlweDimension(1),
        STD_DEV_4096,
        DecompositionBaseLog(11),
        DecompositionLevelCount(3),
        DecompositionBaseLog(2),
        DecompositionLevelCount(8),
        DecompositionBaseLog(12),
        DecompositionLevelCount(4),
        FftType::Split(40),
        DecompositionBaseLog(10),
        DecompositionLevelCount(4),
        DecompositionBaseLog(5),
        DecompositionLevelCount(4),
        LutCountLog(2),
        CiphertextModulus::<u64>::new_native(),
        5,
    );
}

pub struct FheKeys {
    pub lwe_secret_key: LweSecretKey<Vec<u64>>,
    pub glwe_secret_key: GlweSecretKey<Vec<u64>>,
    pub lwe_secret_key_after_ks: LweSecretKey<Vec<u64>>,
    pub bsk: FourierLweBootstrapKey<ABox<[c64]>>,
    pub ksk: LweKeyswitchKey<Vec<u64>>,
    pub auto_keys: HashMap<usize, AutomorphKey<ABox<[c64]>>>,
    pub ss_key: FourierGgswCiphertextList<Vec<c64>>,
}

pub struct FheParams {
    pub log_n: usize,
    pub n: u64,
    pub m: usize,
    pub n_nibbles: usize,
    pub n_rems: usize,
    pub polynomial_size: PolynomialSize,
    pub glwe_size: GlweSize,
    pub glwe_modular_std_dev: StandardDev,
    pub cbs_base_log: DecompositionBaseLog,
    pub cbs_level: DecompositionLevelCount,
    pub ciphertext_modulus: CiphertextModulus<u64>,
    pub encryption_generator: EncryptionRandomGenerator::<ActivatedRandomGenerator>,
    pub delta: u64,
}
