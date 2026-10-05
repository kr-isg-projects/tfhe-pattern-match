use crate::param;
use aligned_vec::ABox;
use dyn_stack::{PodBuffer, PodStack, alloc::Global};
use tfhe::core_crypto::entities::FourierLweBootstrapKey;
use tfhe::core_crypto::fft_impl::fft64::c64;
use tfhe::core_crypto::fft_impl::fft64::crypto::wop_pbs::{
    circuit_bootstrap_boolean, circuit_bootstrap_boolean_scratch, extract_bits,
    extract_bits_scratch,
};
use tfhe::core_crypto::fft_impl::fft64::math::fft::Fft;
use tfhe::core_crypto::prelude::*;

pub fn cbs_int_4to2(
    lwe: &LweCiphertext<Vec<u64>>,
    ksk: &LweKeyswitchKey<Vec<u64>>,
    bsk: &FourierLweBootstrapKey<ABox<[c64]>>,
    pfpksk: &LwePrivateFunctionalPackingKeyswitchKeyList<Vec<u64>>,
    num_ggsw: usize,
) -> GgswCiphertextList<Vec<u64>> {
    let message_size = param::message_size();
    let glwe_dimension = param::glwe_dimension();
    let polynomial_size = param::polynomial_size();
    let cbs_base_log = param::cbs_base_log();
    let cbs_level = param::cbs_level();
    let ciphertext_modulus = param::ciphertext_modulus();
    let glwe_size = glwe_dimension.to_glwe_size();
    let fft = Fft::new(polynomial_size);
    let fft = fft.as_view();
    let delta_log = DeltaLog(u64::BITS as usize - message_size);
    let boolean_delta_log = DeltaLog(u64::BITS as usize - 1);

    let extract_req = extract_bits_scratch::<u64>(
        lwe.lwe_size().to_lwe_dimension(),
        ksk.output_key_lwe_dimension(),
        glwe_size,
        polynomial_size,
        fft,
    );
    let mut extract_mem = PodBuffer::new_in(extract_req, Global);
    let mut extract_stack = PodStack::new(&mut extract_mem);

    let mut extracted_bits = LweCiphertextList::new(
        u64::ZERO,
        ksk.output_lwe_size(),
        LweCiphertextCount(num_ggsw),
        ciphertext_modulus,
    );
    extract_bits(
        extracted_bits.as_mut_view(),
        lwe.as_view(),
        ksk.as_view(),
        bsk.as_view(),
        delta_log,
        ExtractedBitsCount(num_ggsw),
        fft,
        &mut extract_stack,
    );

    let mut ggsw_list_out = GgswCiphertextList::new(
        0u64,
        glwe_size,
        polynomial_size,
        cbs_base_log,
        cbs_level,
        GgswCiphertextCount(num_ggsw),
        ciphertext_modulus,
    );
    let cbs_req = circuit_bootstrap_boolean_scratch::<u64>(
        ksk.output_lwe_size(),
        bsk.output_lwe_dimension().to_lwe_size(),
        glwe_size,
        polynomial_size,
        fft,
    );
    let mut cbs_mem = PodBuffer::new_in(cbs_req, Global);

    for (bit_ct, mut ggsw_out) in extracted_bits.iter().rev().zip(ggsw_list_out.iter_mut()) {
        let mut stack = PodStack::new(&mut cbs_mem);
        circuit_bootstrap_boolean(
            bsk.as_view(),
            bit_ct.as_view(),
            ggsw_out.as_mut_view(),
            boolean_delta_log,
            pfpksk.as_view(),
            fft,
            &mut stack,
        );
    }

    ggsw_list_out
}

pub fn cbs_int_4to1(
    lwe: &LweCiphertext<Vec<u64>>,
    ksk: &LweKeyswitchKey<Vec<u64>>,
    bsk: &FourierLweBootstrapKey<ABox<[c64]>>,
    pfpksk: &LwePrivateFunctionalPackingKeyswitchKeyList<Vec<u64>>,
    idx: usize,
) -> GgswCiphertext<Vec<u64>> {
    let message_size = param::message_size();
    let glwe_dimension = param::glwe_dimension();
    let polynomial_size = param::polynomial_size();
    let cbs_base_log = param::cbs_base_log();
    let cbs_level = param::cbs_level();
    let ciphertext_modulus = param::ciphertext_modulus();
    let glwe_size = glwe_dimension.to_glwe_size();
    let fft = Fft::new(polynomial_size);
    let fft = fft.as_view();
    let delta_log = DeltaLog(u64::BITS as usize - message_size);
    let boolean_delta_log = DeltaLog(u64::BITS as usize - 1);

    let extract_req = extract_bits_scratch::<u64>(
        lwe.lwe_size().to_lwe_dimension(),
        ksk.output_key_lwe_dimension(),
        glwe_size,
        polynomial_size,
        fft,
    );
    let mut extract_mem = PodBuffer::new_in(extract_req, Global);
    let mut extract_stack = PodStack::new(&mut extract_mem);

    let mut extracted_bits = LweCiphertextList::new(
        u64::ZERO,
        ksk.output_lwe_size(),
        LweCiphertextCount(idx + 1),
        ciphertext_modulus,
    );
    extract_bits(
        extracted_bits.as_mut_view(),
        lwe.as_view(),
        ksk.as_view(),
        bsk.as_view(),
        delta_log,
        ExtractedBitsCount(idx + 1),
        fft,
        &mut extract_stack,
    );

    let mut ggsw = GgswCiphertext::new(
        0u64,
        glwe_size,
        polynomial_size,
        cbs_base_log,
        cbs_level,
        ciphertext_modulus,
    );
    let cbs_req = circuit_bootstrap_boolean_scratch::<u64>(
        ksk.output_lwe_size(),
        bsk.output_lwe_dimension().to_lwe_size(),
        glwe_size,
        polynomial_size,
        fft,
    );
    let mut cbs_mem = PodBuffer::new_in(cbs_req, Global);
    let mut stack = PodStack::new(&mut cbs_mem);
    let bit_ct = extracted_bits.get(0);
    circuit_bootstrap_boolean(
        bsk.as_view(),
        bit_ct.as_view(),
        ggsw.as_mut_view(),
        boolean_delta_log,
        pfpksk.as_view(),
        fft,
        &mut stack,
    );

    ggsw
}
