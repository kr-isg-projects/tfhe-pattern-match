use aligned_vec::ABox;
use std::collections::HashMap;

use tfhe::core_crypto::fft_impl::fft64::c64;
use tfhe::core_crypto::prelude::*;

use crate::automorphism::AutomorphKey;

pub const STD_DEV_920: StandardDev = StandardDev(6.843470867416183799e-06);
pub const STD_DEV_2048: StandardDev = StandardDev(3.54e-13);

pub fn lwe_dimension() -> LweDimension {
    LweDimension(920)
}

pub fn lwe_modular_std_dev() -> StandardDev {
    STD_DEV_920
}

pub fn polynomial_size() -> PolynomialSize {
    PolynomialSize(2048)
}

pub fn glwe_dimension() -> GlweDimension {
    GlweDimension(1)
}

pub fn glwe_modular_std_dev() -> StandardDev {
    STD_DEV_2048
}

pub fn pbs_base_log() -> DecompositionBaseLog {
    DecompositionBaseLog(11)
}

pub fn pbs_level() -> DecompositionLevelCount {
    DecompositionLevelCount(3)
}

pub fn ks_base_log() -> DecompositionBaseLog {
    DecompositionBaseLog(2)
}

pub fn ks_level() -> DecompositionLevelCount {
    DecompositionLevelCount(7)
}

pub fn auto_base_log() -> DecompositionBaseLog {
    DecompositionBaseLog(12)
}

pub fn auto_level() -> DecompositionLevelCount {
    DecompositionLevelCount(2)
}

pub fn ss_base_log() -> DecompositionBaseLog {
    DecompositionBaseLog(12)
}

pub fn ss_level() -> DecompositionLevelCount {
    DecompositionLevelCount(2)
}

pub fn cbs_base_log() -> DecompositionBaseLog {
    DecompositionBaseLog(2)
}

pub fn cbs_level() -> DecompositionLevelCount {
    DecompositionLevelCount(7)
}

pub fn ciphertext_modulus() -> CiphertextModulus<u64> {
    CiphertextModulus::<u64>::new_native()
}

pub fn message_size() -> usize {
    4
}

pub struct FheKeys {
    pub lwe_secret_key: LweSecretKey<Vec<u64>>,
    pub glwe_secret_key: GlweSecretKey<Vec<u64>>,
    pub lwe_secret_key_after_ks: LweSecretKey<Vec<u64>>,
    pub bsk: FourierLweBootstrapKey<ABox<[c64]>>,
    pub ksk: LweKeyswitchKey<Vec<u64>>,
    pub pfpksk: LwePrivateFunctionalPackingKeyswitchKeyList<Vec<u64>>,
    pub auto_keys: HashMap<usize, AutomorphKey>,
}

pub struct FheParams {
    pub log_n: usize,
    pub n: u64,
    pub m: usize,
    pub n_parts: usize,
    pub top_part_bits: usize,
    pub polynomial_size: PolynomialSize,
    pub glwe_size: GlweSize,
    pub glwe_modular_std_dev: StandardDev,
    pub cbs_base_log: DecompositionBaseLog,
    pub cbs_level: DecompositionLevelCount,
    pub ciphertext_modulus: CiphertextModulus<u64>,
    pub encryption_generator: EncryptionRandomGenerator<DefaultRandomGenerator>,
    pub delta: u64,
}
