use crate::param;
use crate::wop_pbs::{
    private_vertical_packing, private_vertical_packing_multi_extract,
    private_vertical_packing_to_glwe, vertical_packing, vertical_packing_multi_extract,
    vertical_packing_scratch,
};
use tfhe::core_crypto::fft_impl::fft64::c64;
use tfhe::core_crypto::prelude::*;

pub fn lut_eval_n_to_4(
    lut: &PolynomialList<Vec<u64>>, // polynomial is used instead of RLWE
    ggsw_list: &GgswCiphertextList<Vec<u64>>,
    lwe_secret_key: &LweSecretKey<Vec<u64>>,
    polynomial_size: PolynomialSize,
    cbs_base_log: DecompositionBaseLog,
    cbs_level: DecompositionLevelCount,
    ciphertext_modulus: CiphertextModulus<u64>,
) -> LweCiphertext<Vec<u64>> {
    let glwe_dimension = param::glwe_dimension();
    let glwe_size = glwe_dimension.to_glwe_size();
    let lut_input_size = ggsw_list.iter().len();
    let num_lut_in_vp = lut.polynomial_count().0 as usize;
    if lut_input_size < polynomial_size.0.trailing_zeros() as usize {
        unimplemented!(
            "LUT evaluation using only blind rotation for log_n < 12 is not implemented yet"
        );
    }
    let mut fourier_ggsw_list = FourierGgswCiphertextList::new(
        vec![
            c64::default();
            lut_input_size
                * polynomial_size.to_fourier_polynomial_size().0
                * glwe_size.0
                * glwe_size.0
                * cbs_level.0
        ],
        lut_input_size,
        glwe_size,
        polynomial_size,
        cbs_base_log,
        cbs_level,
    );
    for (mut fourier_ggsw, ggsw) in fourier_ggsw_list
        .as_mut_view()
        .into_ggsw_iter()
        .zip(ggsw_list.iter())
    {
        convert_standard_ggsw_ciphertext_to_fourier(&ggsw, &mut fourier_ggsw);
    }
    let fft = Fft::new(polynomial_size);
    let fft = fft.as_view();
    let mut lwe_out = LweCiphertext::new(
        0u64,
        lwe_secret_key.lwe_dimension().to_lwe_size(),
        ciphertext_modulus,
    );
    let mut buffer = ComputationBuffers::new();
    let scratch = vertical_packing_scratch::<u64>(
        glwe_size,
        polynomial_size,
        PolynomialCount(num_lut_in_vp),
        lut_input_size,
        fft,
    );
    buffer.resize(scratch.unaligned_bytes_required());
    let stack = buffer.stack();
    vertical_packing(
        lut.as_view(),
        lwe_out.as_mut_view(),
        fourier_ggsw_list.as_view(),
        fft,
        stack,
    );
    lwe_out
}

pub fn private_lut_eval_n_to_4(
    lut: &GlweCiphertextList<Vec<u64>>,
    ggsw_list: &GgswCiphertextList<Vec<u64>>,
    lwe_secret_key: &LweSecretKey<Vec<u64>>,
    polynomial_size: PolynomialSize,
    cbs_base_log: DecompositionBaseLog,
    cbs_level: DecompositionLevelCount,
    ciphertext_modulus: CiphertextModulus<u64>,
) -> LweCiphertext<Vec<u64>> {
    let glwe_dimension = param::glwe_dimension();
    let glwe_size = glwe_dimension.to_glwe_size();
    let lut_input_size = ggsw_list.iter().len();
    let num_lut_in_vp = lut.glwe_ciphertext_count().0 as usize;
    if lut_input_size < polynomial_size.0.trailing_zeros() as usize {
        unimplemented!(
            "LUT evaluation using only blind rotation for log_n < 12 is not implemented yet"
        );
    }
    let lut_view = GlweCiphertextList::from_container(
        lut.as_ref(),
        lut.glwe_size(),
        lut.polynomial_size(),
        lut.ciphertext_modulus(),
    );
    let mut fourier_ggsw_list = FourierGgswCiphertextList::new(
        vec![
            c64::default();
            lut_input_size
                * polynomial_size.to_fourier_polynomial_size().0
                * glwe_size.0
                * glwe_size.0
                * cbs_level.0
        ],
        lut_input_size,
        glwe_size,
        polynomial_size,
        cbs_base_log,
        cbs_level,
    );
    for (mut fourier_ggsw, ggsw) in fourier_ggsw_list
        .as_mut_view()
        .into_ggsw_iter()
        .zip(ggsw_list.iter())
    {
        convert_standard_ggsw_ciphertext_to_fourier(&ggsw, &mut fourier_ggsw);
    }
    let fft = Fft::new(polynomial_size);
    let fft = fft.as_view();
    let mut lwe_out = LweCiphertext::new(
        0u64,
        lwe_secret_key.lwe_dimension().to_lwe_size(),
        ciphertext_modulus,
    );
    let mut buffer = ComputationBuffers::new();
    let scratch = vertical_packing_scratch::<u64>(
        glwe_size,
        polynomial_size,
        PolynomialCount(num_lut_in_vp),
        lut_input_size,
        fft,
    );
    buffer.resize(scratch.unaligned_bytes_required());
    let stack = buffer.stack();
    private_vertical_packing(
        lut_view,
        lwe_out.as_mut_view(),
        fourier_ggsw_list.as_view(),
        fft,
        stack,
    );
    lwe_out
}

pub fn lut_eval_n_to_4_slice_multi(
    lut: &PolynomialList<Vec<u64>>, // polynomial is used instead of RLWE
    ggsw_list: &GgswCiphertextList<Vec<u64>>,
    lwe_secret_key: &LweSecretKey<Vec<u64>>,
    polynomial_size: PolynomialSize,
    cbs_base_log: DecompositionBaseLog,
    cbs_level: DecompositionLevelCount,
    ciphertext_modulus: CiphertextModulus<u64>,
    extract_len: usize,
) -> Vec<LweCiphertext<Vec<u64>>> {
    let glwe_dimension = param::glwe_dimension();
    let glwe_size = glwe_dimension.to_glwe_size();
    let lut_input_size = ggsw_list.iter().len();
    let num_lut_in_vp = lut.polynomial_count().0 as usize;
    if lut_input_size < polynomial_size.0.trailing_zeros() as usize {
        unimplemented!(
            "LUT evaluation using only blind rotation for log_n < 12 is not implemented yet"
        );
    }
    let mut fourier_ggsw_list = FourierGgswCiphertextList::new(
        vec![
            c64::default();
            lut_input_size
                * polynomial_size.to_fourier_polynomial_size().0
                * glwe_size.0
                * glwe_size.0
                * cbs_level.0
        ],
        lut_input_size,
        glwe_size,
        polynomial_size,
        cbs_base_log,
        cbs_level,
    );
    for (mut fourier_ggsw, ggsw) in fourier_ggsw_list
        .as_mut_view()
        .into_ggsw_iter()
        .zip(ggsw_list.iter())
    {
        convert_standard_ggsw_ciphertext_to_fourier(&ggsw, &mut fourier_ggsw);
    }
    let fft = Fft::new(polynomial_size);
    let fft = fft.as_view();
    let mut buffer = ComputationBuffers::new();
    let scratch = vertical_packing_scratch::<u64>(
        glwe_size,
        polynomial_size,
        PolynomialCount(num_lut_in_vp),
        lut_input_size,
        fft,
    );
    buffer.resize(scratch.unaligned_bytes_required());
    let stack = buffer.stack();
    let mut lwe_outputs = Vec::with_capacity(extract_len);
    for _ in 0..extract_len {
        lwe_outputs.push(LweCiphertext::new(
            0u64,
            lwe_secret_key.lwe_dimension().to_lwe_size(),
            ciphertext_modulus,
        ));
    }
    vertical_packing_multi_extract(
        lut.as_view(),
        fourier_ggsw_list.as_view(),
        fft,
        stack,
        &mut lwe_outputs,
        ciphertext_modulus,
        lwe_secret_key.lwe_dimension(),
    );
    lwe_outputs
}

pub fn private_lut_eval_n_to_4_slice_multi(
    lut: &GlweCiphertextList<Vec<u64>>,
    ggsw_list: &GgswCiphertextList<Vec<u64>>,
    lwe_secret_key: &LweSecretKey<Vec<u64>>,
    polynomial_size: PolynomialSize,
    cbs_base_log: DecompositionBaseLog,
    cbs_level: DecompositionLevelCount,
    ciphertext_modulus: CiphertextModulus<u64>,
    extract_len: usize,
) -> Vec<LweCiphertext<Vec<u64>>> {
    let glwe_dimension = param::glwe_dimension();
    let glwe_size = glwe_dimension.to_glwe_size();
    let lut_input_size = ggsw_list.iter().len();
    let num_lut_in_vp = lut.glwe_ciphertext_count().0 as usize;
    if lut_input_size < polynomial_size.0.trailing_zeros() as usize {
        unimplemented!(
            "LUT evaluation using only blind rotation for log_n < 12 is not implemented yet"
        );
    }
    let lut_view = GlweCiphertextList::from_container(
        lut.as_ref(),
        lut.glwe_size(),
        lut.polynomial_size(),
        lut.ciphertext_modulus(),
    );
    let mut fourier_ggsw_list = FourierGgswCiphertextList::new(
        vec![
            c64::default();
            lut_input_size
                * polynomial_size.to_fourier_polynomial_size().0
                * glwe_size.0
                * glwe_size.0
                * cbs_level.0
        ],
        lut_input_size,
        glwe_size,
        polynomial_size,
        cbs_base_log,
        cbs_level,
    );
    for (mut fourier_ggsw, ggsw) in fourier_ggsw_list
        .as_mut_view()
        .into_ggsw_iter()
        .zip(ggsw_list.iter())
    {
        convert_standard_ggsw_ciphertext_to_fourier(&ggsw, &mut fourier_ggsw);
    }
    let fft = Fft::new(polynomial_size);
    let fft = fft.as_view();
    let mut buffer = ComputationBuffers::new();
    let scratch = vertical_packing_scratch::<u64>(
        glwe_size,
        polynomial_size,
        PolynomialCount(num_lut_in_vp),
        lut_input_size,
        fft,
    );
    buffer.resize(scratch.unaligned_bytes_required());
    let stack = buffer.stack();
    let mut lwe_outputs = Vec::with_capacity(extract_len);
    for _ in 0..extract_len {
        lwe_outputs.push(LweCiphertext::new(
            0u64,
            lwe_secret_key.lwe_dimension().to_lwe_size(),
            ciphertext_modulus,
        ));
    }
    private_vertical_packing_multi_extract(
        lut_view,
        fourier_ggsw_list.as_view(),
        fft,
        stack,
        &mut lwe_outputs,
        ciphertext_modulus,
        lwe_secret_key.lwe_dimension(),
    );
    lwe_outputs
}

pub fn private_lut_eval_n_to_4_rotated_rlwe(
    lut: &GlweCiphertextList<Vec<u64>>,
    ggsw_list: &GgswCiphertextList<Vec<u64>>,
    polynomial_size: PolynomialSize,
    cbs_base_log: DecompositionBaseLog,
    cbs_level: DecompositionLevelCount,
    ciphertext_modulus: CiphertextModulus<u64>,
) -> GlweCiphertext<Vec<u64>> {
    let glwe_dimension = param::glwe_dimension();
    let glwe_size = glwe_dimension.to_glwe_size();
    let lut_input_size = ggsw_list.iter().len();
    let num_lut_in_vp = lut.glwe_ciphertext_count().0 as usize;
    if lut_input_size < polynomial_size.0.trailing_zeros() as usize {
        unimplemented!(
            "LUT evaluation using only blind rotation for log_n < 12 is not implemented yet"
        );
    }
    let lut_view = GlweCiphertextList::from_container(
        lut.as_ref(),
        lut.glwe_size(),
        lut.polynomial_size(),
        lut.ciphertext_modulus(),
    );
    let mut fourier_ggsw_list = FourierGgswCiphertextList::new(
        vec![
            c64::default();
            lut_input_size
                * polynomial_size.to_fourier_polynomial_size().0
                * glwe_size.0
                * glwe_size.0
                * cbs_level.0
        ],
        lut_input_size,
        glwe_size,
        polynomial_size,
        cbs_base_log,
        cbs_level,
    );
    for (mut fourier_ggsw, ggsw) in fourier_ggsw_list
        .as_mut_view()
        .into_ggsw_iter()
        .zip(ggsw_list.iter())
    {
        convert_standard_ggsw_ciphertext_to_fourier(&ggsw, &mut fourier_ggsw);
    }
    let fft = Fft::new(polynomial_size);
    let fft = fft.as_view();
    let mut buffer = ComputationBuffers::new();
    let scratch = vertical_packing_scratch::<u64>(
        glwe_size,
        polynomial_size,
        PolynomialCount(num_lut_in_vp),
        lut_input_size,
        fft,
    );
    buffer.resize(scratch.unaligned_bytes_required());
    let stack = buffer.stack();

    let mut glwe_out = GlweCiphertext::new(0u64, glwe_size, polynomial_size, ciphertext_modulus);
    private_vertical_packing_to_glwe(
        lut_view,
        glwe_out.as_mut_view(),
        fourier_ggsw_list.as_view(),
        fft,
        stack,
    );
    glwe_out
}
