use dyn_stack::PodStack;

use tfhe::core_crypto::fft_impl::fft64::crypto::ggsw::FourierGgswCiphertextListView;
use tfhe::core_crypto::fft_impl::fft64::crypto::ggsw::add_external_product_assign;
use tfhe::core_crypto::fft_impl::fft64::crypto::wop_pbs::{
    blind_rotate_assign, cmux_tree_memory_optimized, cmux_tree_memory_optimized_scratch,
};
pub use tfhe::core_crypto::fft_impl::fft64::crypto::wop_pbs::{
    vertical_packing, vertical_packing_scratch,
};
use tfhe::core_crypto::fft_impl::fft64::math::fft::FftView;
use tfhe::core_crypto::prelude::*;

pub fn vertical_packing_multi_extract<Scalar: UnsignedTorus + CastInto<usize>>(
    lut: PolynomialList<&[Scalar]>,
    ggsw_list: FourierGgswCiphertextListView<'_>,
    fft: FftView<'_>,
    stack: &mut PodStack,
    lwe_outputs: &mut [LweCiphertext<Vec<Scalar>>],
    ciphertext_modulus: CiphertextModulus<Scalar>,
    _lwe_dimension: LweDimension,
) {
    let polynomial_size = ggsw_list.polynomial_size();
    let glwe_size = ggsw_list.glwe_size();

    let log_lut_number: usize =
        Scalar::BITS - 1 - lut.polynomial_count().0.leading_zeros() as usize;

    let log_number_of_luts_for_cmux_tree = if log_lut_number > ggsw_list.count() {
        0
    } else {
        log_lut_number
    };

    let (cmux_ggsw, br_ggsw) = ggsw_list.split_at(log_number_of_luts_for_cmux_tree);

    let (cmux_tree_lut_res_data, stack) = stack.make_aligned_with(
        polynomial_size.0 * glwe_size.0,
        aligned_vec::CACHELINE_ALIGN,
        |_| Scalar::ZERO,
    );
    let mut cmux_tree_lut_res = GlweCiphertext::from_container(
        &mut *cmux_tree_lut_res_data,
        polynomial_size,
        ciphertext_modulus,
    );

    cmux_tree_memory_optimized(cmux_tree_lut_res.as_mut_view(), lut, cmux_ggsw, fft, stack);
    blind_rotate_assign(cmux_tree_lut_res.as_mut_view(), br_ggsw, fft, stack);

    for (i, lwe_out) in lwe_outputs.iter_mut().enumerate() {
        extract_lwe_sample_from_glwe_ciphertext(&cmux_tree_lut_res, lwe_out, MonomialDegree(i));
    }
}

pub fn private_cmux_tree<Scalar: UnsignedTorus + CastInto<usize>>(
    mut output_glwe: GlweCiphertext<&mut [Scalar]>,
    lut_per_layer: GlweCiphertextList<&[Scalar]>,
    ggsw_list: FourierGgswCiphertextListView<'_>,
    fft: FftView<'_>,
    stack: &mut PodStack,
) {
    if ggsw_list.count() > 0 {
        let glwe_size = output_glwe.glwe_size();
        let ciphertext_modulus = output_glwe.ciphertext_modulus();
        let polynomial_size = ggsw_list.polynomial_size();
        let nb_layer = ggsw_list.count();

        debug_assert!(stack.can_hold(cmux_tree_memory_optimized_scratch::<Scalar>(
            glwe_size,
            polynomial_size,
            nb_layer,
            fft
        )));

        let (t_0_data, stack) = stack.make_aligned_with(
            polynomial_size.0 * glwe_size.0 * nb_layer,
            aligned_vec::CACHELINE_ALIGN,
            |_| Scalar::ZERO,
        );
        let (t_1_data, stack) = stack.make_aligned_with(
            polynomial_size.0 * glwe_size.0 * nb_layer,
            aligned_vec::CACHELINE_ALIGN,
            |_| Scalar::ZERO,
        );

        let mut t_0 = GlweCiphertextList::from_container(
            t_0_data.as_mut(),
            glwe_size,
            polynomial_size,
            ciphertext_modulus,
        );
        let mut t_1 = GlweCiphertextList::from_container(
            t_1_data.as_mut(),
            glwe_size,
            polynomial_size,
            ciphertext_modulus,
        );

        let (t_fill, stack) = stack.make_with(nb_layer, |_| 0_usize);

        let mut lut_glwe_iter = lut_per_layer.iter();
        loop {
            let even = lut_glwe_iter.next();
            let odd = lut_glwe_iter.next();

            let (Some(lut_2i), Some(lut_2i_plus_1)) = (even, odd) else {
                break;
            };

            let mut t_iter = t_0.iter_mut().zip(t_1.iter_mut()).enumerate();
            let (mut j_counter, (mut t0_j, mut t1_j)) = t_iter.next().unwrap();

            t0_j.as_mut().copy_from_slice(lut_2i.as_ref());
            t1_j.as_mut().copy_from_slice(lut_2i_plus_1.as_ref());

            t_fill[0] = 2;

            for (j, ggsw) in ggsw_list.into_ggsw_iter().rev().enumerate() {
                if t_fill[j] == 2 {
                    let (diff_data, mut stack_local) = stack.collect_aligned(
                        aligned_vec::CACHELINE_ALIGN,
                        t1_j.as_ref()
                            .iter()
                            .zip(t0_j.as_ref().iter())
                            .map(|(&a, &b)| a.wrapping_sub(b)),
                    );
                    let diff = GlweCiphertext::from_container(
                        &*diff_data,
                        polynomial_size,
                        ciphertext_modulus,
                    );

                    if j < nb_layer - 1 {
                        let (j_counter_plus_1, (mut t_0_j_plus_1, mut t_1_j_plus_1)) =
                            t_iter.next().unwrap();

                        assert_eq!(j_counter, j);
                        assert_eq!(j_counter_plus_1, j + 1);

                        let mut output = if t_fill[j + 1] == 0 {
                            t_0_j_plus_1.as_mut_view()
                        } else {
                            t_1_j_plus_1.as_mut_view()
                        };

                        output.as_mut().copy_from_slice(t0_j.as_ref());
                        add_external_product_assign(output, ggsw, diff, fft, &mut stack_local);
                        t_fill[j + 1] += 1;
                        t_fill[j] = 0;

                        (j_counter, t0_j, t1_j) = (j_counter_plus_1, t_0_j_plus_1, t_1_j_plus_1);
                    } else {
                        assert_eq!(j, nb_layer - 1);
                        let mut output = output_glwe.as_mut_view();
                        output.as_mut().copy_from_slice(t0_j.as_ref());
                        add_external_product_assign(output, ggsw, diff, fft, &mut stack_local);
                    }
                } else {
                    break;
                }
            }
        }
    } else {
        output_glwe.as_mut().copy_from_slice(lut_per_layer.as_ref());
    }
}

pub fn private_vertical_packing<Scalar: UnsignedTorus + CastInto<usize>>(
    lut: GlweCiphertextList<&[Scalar]>,
    mut lwe_out: LweCiphertext<&mut [Scalar]>,
    ggsw_list: FourierGgswCiphertextListView<'_>,
    fft: FftView<'_>,
    stack: &mut PodStack,
) {
    debug_assert!(
        lwe_out.ciphertext_modulus().is_native_modulus(),
        "This operation currently only supports native moduli"
    );

    let polynomial_size = ggsw_list.polynomial_size();
    let glwe_size = ggsw_list.glwe_size();
    let glwe_dimension = glwe_size.to_glwe_dimension();
    let ciphertext_modulus = lwe_out.ciphertext_modulus();

    debug_assert!(
        lwe_out.lwe_size().to_lwe_dimension()
            == glwe_dimension.to_equivalent_lwe_dimension(polynomial_size),
        "Output LWE ciphertext needs to have an LweDimension of {:?}, got {:?}",
        glwe_dimension.to_equivalent_lwe_dimension(polynomial_size),
        lwe_out.lwe_size().to_lwe_dimension(),
    );

    let log_lut_number: usize =
        Scalar::BITS - 1 - lut.glwe_ciphertext_count().0.leading_zeros() as usize;

    let log_number_of_luts_for_cmux_tree = if log_lut_number > ggsw_list.count() {
        0
    } else {
        log_lut_number
    };

    let (cmux_ggsw, br_ggsw) = ggsw_list.split_at(log_number_of_luts_for_cmux_tree);

    let (cmux_tree_lut_res_data, stack) = stack.make_aligned_with(
        polynomial_size.0 * glwe_size.0,
        aligned_vec::CACHELINE_ALIGN,
        |_| Scalar::ZERO,
    );
    let mut cmux_tree_lut_res = GlweCiphertext::from_container(
        &mut *cmux_tree_lut_res_data,
        polynomial_size,
        ciphertext_modulus,
    );

    private_cmux_tree(cmux_tree_lut_res.as_mut_view(), lut, cmux_ggsw, fft, stack);
    blind_rotate_assign(cmux_tree_lut_res.as_mut_view(), br_ggsw, fft, stack);

    extract_lwe_sample_from_glwe_ciphertext(&cmux_tree_lut_res, &mut lwe_out, MonomialDegree(0));
}

pub fn private_vertical_packing_multi_extract<Scalar: UnsignedTorus + CastInto<usize>>(
    lut: GlweCiphertextList<&[Scalar]>,
    ggsw_list: FourierGgswCiphertextListView<'_>,
    fft: FftView<'_>,
    stack: &mut PodStack,
    lwe_outputs: &mut [LweCiphertext<Vec<Scalar>>],
    ciphertext_modulus: CiphertextModulus<Scalar>,
    _lwe_dimension: LweDimension,
) {
    let polynomial_size = ggsw_list.polynomial_size();
    let glwe_size = ggsw_list.glwe_size();

    let log_lut_number: usize =
        Scalar::BITS - 1 - lut.glwe_ciphertext_count().0.leading_zeros() as usize;

    let log_number_of_luts_for_cmux_tree = if log_lut_number > ggsw_list.count() {
        0
    } else {
        log_lut_number
    };

    let (cmux_ggsw, br_ggsw) = ggsw_list.split_at(log_number_of_luts_for_cmux_tree);

    let (cmux_tree_lut_res_data, stack) = stack.make_aligned_with(
        polynomial_size.0 * glwe_size.0,
        aligned_vec::CACHELINE_ALIGN,
        |_| Scalar::ZERO,
    );
    let mut cmux_tree_lut_res = GlweCiphertext::from_container(
        &mut *cmux_tree_lut_res_data,
        polynomial_size,
        ciphertext_modulus,
    );

    private_cmux_tree(cmux_tree_lut_res.as_mut_view(), lut, cmux_ggsw, fft, stack);
    blind_rotate_assign(cmux_tree_lut_res.as_mut_view(), br_ggsw, fft, stack);

    for (i, lwe_out) in lwe_outputs.iter_mut().enumerate() {
        extract_lwe_sample_from_glwe_ciphertext(&cmux_tree_lut_res, lwe_out, MonomialDegree(i));
    }
}

pub fn private_vertical_packing_to_glwe<Scalar: UnsignedTorus + CastInto<usize>>(
    lut: GlweCiphertextList<&[Scalar]>,
    mut glwe_out: GlweCiphertext<&mut [Scalar]>,
    ggsw_list: FourierGgswCiphertextListView<'_>,
    fft: FftView<'_>,
    stack: &mut PodStack,
) {
    let polynomial_size = ggsw_list.polynomial_size();
    let glwe_size = ggsw_list.glwe_size();
    let ciphertext_modulus = glwe_out.ciphertext_modulus();

    debug_assert_eq!(glwe_out.glwe_size(), glwe_size);
    debug_assert_eq!(glwe_out.polynomial_size(), polynomial_size);

    let log_lut_number: usize =
        Scalar::BITS - 1 - lut.glwe_ciphertext_count().0.leading_zeros() as usize;

    let log_number_of_luts_for_cmux_tree = if log_lut_number > ggsw_list.count() {
        0
    } else {
        log_lut_number
    };

    let (cmux_ggsw, br_ggsw) = ggsw_list.split_at(log_number_of_luts_for_cmux_tree);

    let (cmux_tree_lut_res_data, stack) = stack.make_aligned_with(
        polynomial_size.0 * glwe_size.0,
        aligned_vec::CACHELINE_ALIGN,
        |_| Scalar::ZERO,
    );
    let mut cmux_tree_lut_res = GlweCiphertext::from_container(
        &mut *cmux_tree_lut_res_data,
        polynomial_size,
        ciphertext_modulus,
    );

    private_cmux_tree(cmux_tree_lut_res.as_mut_view(), lut, cmux_ggsw, fft, stack);
    blind_rotate_assign(cmux_tree_lut_res.as_mut_view(), br_ggsw, fft, stack);

    glwe_out
        .as_mut()
        .copy_from_slice(cmux_tree_lut_res.as_ref());
}
