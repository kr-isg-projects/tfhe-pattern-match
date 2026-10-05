use std::collections::HashMap;

use tfhe::core_crypto::{algorithms::slice_algorithms::slice_wrapping_opposite_assign, prelude::*};

use crate::automorphism::AutomorphKey;

fn glwe_ciphertext_clone_from<OutputCont, InputCont>(
    dst: &mut GlweCiphertext<OutputCont>,
    src: &GlweCiphertext<InputCont>,
) where
    InputCont: Container<Element = u64>,
    OutputCont: ContainerMut<Element = u64>,
{
    dst.as_mut().clone_from_slice(src.as_ref());
}

fn lwe_ciphertext_mod_switch_from_native_to_non_native_power_of_two<InputCont, OutputCont>(
    input: &LweCiphertext<InputCont>,
    output: &mut LweCiphertext<OutputCont>,
) where
    InputCont: Container<Element = u64>,
    OutputCont: ContainerMut<Element = u64>,
{
    let divisor = output
        .ciphertext_modulus()
        .get_power_of_two_scaling_to_native_torus();
    for (src, dst) in input.as_ref().iter().zip(output.as_mut().iter_mut()) {
        *dst = *src - (*src) % divisor;
    }
}

fn lwe_ciphertext_mod_raise_from_non_native_power_of_two_to_native<InputCont, OutputCont>(
    input: &LweCiphertext<InputCont>,
    output: &mut LweCiphertext<OutputCont>,
) where
    InputCont: Container<Element = u64>,
    OutputCont: ContainerMut<Element = u64>,
{
    let scaling_factor = input
        .ciphertext_modulus()
        .get_power_of_two_scaling_to_native_torus();
    for (src, dst) in input.as_ref().iter().zip(output.as_mut().iter_mut()) {
        *dst = src.wrapping_div(scaling_factor);
    }
}

fn lwe_preprocessing_assign<ContMut>(
    input: &mut LweCiphertext<ContMut>,
    polynomial_size: PolynomialSize,
) where
    ContMut: ContainerMut<Element = u64>,
{
    let log_small_q = u64::BITS as usize - polynomial_size.0.ilog2() as usize;
    let small_ciphertext_modulus =
        CiphertextModulus::<u64>::try_new_power_of_2(log_small_q).unwrap();
    let mut buf = LweCiphertext::new(0u64, input.lwe_size(), small_ciphertext_modulus);

    lwe_ciphertext_mod_switch_from_native_to_non_native_power_of_two(input, &mut buf);
    lwe_ciphertext_mod_raise_from_non_native_power_of_two_to_native(&buf, input);
}

fn lwe_preprocessing<InputCont, OutputCont>(
    input: &LweCiphertext<InputCont>,
    output: &mut LweCiphertext<OutputCont>,
    polynomial_size: PolynomialSize,
) where
    InputCont: Container<Element = u64>,
    OutputCont: ContainerMut<Element = u64>,
{
    output.as_mut().clone_from_slice(input.as_ref());
    lwe_preprocessing_assign(output, polynomial_size);
}

fn convert_lwe_to_glwe_const<InputCont, OutputCont>(
    input: &LweCiphertext<InputCont>,
    output: &mut GlweCiphertext<OutputCont>,
) where
    InputCont: Container<Element = u64>,
    OutputCont: ContainerMut<Element = u64>,
{
    let lwe_dimension = input.lwe_size().to_lwe_dimension().0;
    let glwe_dimension = output.glwe_size().to_glwe_dimension().0;
    let polynomial_size = output.polynomial_size().0;
    let (lwe_mask, lwe_body) = input.get_mask_and_body();
    let (mut glwe_mask, mut glwe_body) = output.get_mut_mask_and_body();

    assert_eq!(lwe_dimension, glwe_dimension * polynomial_size);

    *glwe_body.as_mut().get_mut(0).unwrap() = *lwe_body.data;

    for (glwe_poly, lwe_poly) in glwe_mask
        .as_mut()
        .chunks_exact_mut(polynomial_size)
        .zip(lwe_mask.as_ref().chunks_exact(polynomial_size))
    {
        glwe_poly.clone_from_slice(lwe_poly);
        glwe_poly.reverse();
        slice_wrapping_opposite_assign(&mut glwe_poly[0..(polynomial_size - 1)]);
        glwe_poly.rotate_left(polynomial_size - 1);
    }
}

pub fn convert_lwe_to_glwe_by_trace_with_preprocessing<InputCont, OutputCont>(
    input: &LweCiphertext<InputCont>,
    output: &mut GlweCiphertext<OutputCont>,
    auto_keys: &HashMap<usize, AutomorphKey>,
) where
    InputCont: Container<Element = u64>,
    OutputCont: ContainerMut<Element = u64>,
{
    let lwe_size = input.lwe_size();
    let polynomial_size = output.polynomial_size();
    let ciphertext_modulus = input.ciphertext_modulus();

    let mut buf = LweCiphertext::new(0u64, lwe_size, ciphertext_modulus);
    lwe_preprocessing(input, &mut buf, polynomial_size);
    convert_lwe_to_glwe_const(&buf, output);
    trace_assign(output, auto_keys);
}

fn trace_assign<ContMut>(
    glwe_in: &mut GlweCiphertext<ContMut>,
    auto_keys: &HashMap<usize, AutomorphKey>,
) where
    ContMut: ContainerMut<Element = u64>,
{
    trace_partial_assign(glwe_in, auto_keys, 1);
}

fn trace_partial_assign<ContMut>(
    input: &mut GlweCiphertext<ContMut>,
    auto_keys: &HashMap<usize, AutomorphKey>,
    n: usize,
) where
    ContMut: ContainerMut<Element = u64>,
{
    let glwe_size = input.glwe_size();
    let polynomial_size = input.polynomial_size();
    let ciphertext_modulus = input.ciphertext_modulus();

    let mut buf = GlweCiphertextOwned::new(0u64, glwe_size, polynomial_size, ciphertext_modulus);
    let mut out = GlweCiphertext::new(0u64, glwe_size, polynomial_size, ciphertext_modulus);
    glwe_ciphertext_clone_from(&mut out, input);

    let log_polynomial_size = polynomial_size.0.ilog2() as usize;
    let log_n = n.ilog2() as usize;
    for i in 1..=(log_polynomial_size - log_n) {
        let k = polynomial_size.0 / (1 << (i - 1)) + 1;
        let auto_key = auto_keys.get(&k).unwrap();
        auto_key.auto(&mut buf, &out);
        glwe_ciphertext_add_assign(&mut out, &buf);
    }

    glwe_ciphertext_clone_from(input, &out);
}
