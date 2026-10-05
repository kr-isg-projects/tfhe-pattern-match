use std::collections::HashMap;

use tfhe::core_crypto::{
    algorithms::slice_algorithms::slice_wrapping_opposite_assign,
    prelude::{polynomial_algorithms::*, *},
};

fn eval_x_k(poly: PolynomialView<'_, u64>, k: usize) -> PolynomialOwned<u64> {
    let mut out = PolynomialOwned::new(0u64, poly.polynomial_size());
    out.as_mut()[0] = poly.as_ref()[0];
    for i in 1..poly.polynomial_size().0 {
        let j = i * k % poly.polynomial_size().0;
        let sign = if ((i * k) / poly.polynomial_size().0).is_multiple_of(2) {
            1u64
        } else {
            u64::MAX
        };
        out.as_mut()[j] = sign.wrapping_mul(poly.as_ref()[i]);
    }
    out
}

#[derive(Clone)]
pub struct GlweKeyswitchKey {
    data: Vec<u64>,
    output_glwe_dimension: GlweDimension,
    polynomial_size: PolynomialSize,
    decomp_base_log: DecompositionBaseLog,
    decomp_level_count: DecompositionLevelCount,
    ciphertext_modulus: CiphertextModulus<u64>,
}

impl GlweKeyswitchKey {
    fn new(
        input_glwe_dimension: GlweDimension,
        output_glwe_dimension: GlweDimension,
        polynomial_size: PolynomialSize,
        decomp_base_log: DecompositionBaseLog,
        decomp_level_count: DecompositionLevelCount,
        ciphertext_modulus: CiphertextModulus<u64>,
    ) -> Self {
        let output_glwe_size = output_glwe_dimension.to_glwe_size();
        Self {
            data: vec![
                0u64;
                input_glwe_dimension.0
                    * output_glwe_size.0
                    * polynomial_size.0
                    * decomp_level_count.0
            ],
            output_glwe_dimension,
            polynomial_size,
            decomp_base_log,
            decomp_level_count,
            ciphertext_modulus,
        }
    }

    fn output_glwe_dimension(&self) -> GlweDimension {
        self.output_glwe_dimension
    }

    fn polynomial_size(&self) -> PolynomialSize {
        self.polynomial_size
    }

    fn decomp_base_log(&self) -> DecompositionBaseLog {
        self.decomp_base_log
    }

    fn decomp_level_count(&self) -> DecompositionLevelCount {
        self.decomp_level_count
    }

    fn ciphertext_modulus(&self) -> CiphertextModulus<u64> {
        self.ciphertext_modulus
    }

    fn output_glwe_size(&self) -> GlweSize {
        self.output_glwe_dimension.to_glwe_size()
    }

    fn glev_chunk_len(&self) -> usize {
        self.output_glwe_size().0 * self.polynomial_size.0 * self.decomp_level_count.0
    }

    fn glev_list(&self, index: usize) -> GlweCiphertextList<&[u64]> {
        let chunk_len = self.glev_chunk_len();
        let start = index * chunk_len;
        let end = start + chunk_len;
        GlweCiphertextList::from_container(
            &self.data[start..end],
            self.output_glwe_size(),
            self.polynomial_size,
            self.ciphertext_modulus,
        )
    }

    fn glev_list_mut(&mut self, index: usize) -> GlweCiphertextList<&mut [u64]> {
        let chunk_len = self.glev_chunk_len();
        let start = index * chunk_len;
        let end = start + chunk_len;
        let output_glwe_size = self.output_glwe_size();
        let polynomial_size = self.polynomial_size;
        let ciphertext_modulus = self.ciphertext_modulus;
        GlweCiphertextList::from_container(
            &mut self.data[start..end],
            output_glwe_size,
            polynomial_size,
            ciphertext_modulus,
        )
    }
}

fn allocate_and_generate_new_glwe_keyswitch_key<G: ByteRandomGenerator>(
    input_glwe_sk: &GlweSecretKey<Vec<u64>>,
    output_glwe_sk: &GlweSecretKey<Vec<u64>>,
    decomp_base_log: DecompositionBaseLog,
    decomp_level_count: DecompositionLevelCount,
    noise_parameters: impl DispersionParameter,
    ciphertext_modulus: CiphertextModulus<u64>,
    generator: &mut EncryptionRandomGenerator<G>,
) -> GlweKeyswitchKey {
    let mut glwe_keyswitch_key = GlweKeyswitchKey::new(
        input_glwe_sk.glwe_dimension(),
        output_glwe_sk.glwe_dimension(),
        input_glwe_sk.polynomial_size(),
        decomp_base_log,
        decomp_level_count,
        ciphertext_modulus,
    );

    for (index, input_sk_poly) in input_glwe_sk.as_polynomial_list().iter().enumerate() {
        let mut neg_sk_poly =
            PlaintextList::new(0u64, PlaintextCount(input_sk_poly.polynomial_size().0));
        neg_sk_poly
            .as_mut()
            .clone_from_slice(input_sk_poly.as_ref());
        slice_wrapping_opposite_assign(neg_sk_poly.as_mut());

        for (level, mut glwe) in glwe_keyswitch_key
            .glev_list_mut(index)
            .iter_mut()
            .enumerate()
        {
            let log_scale = u64::BITS as usize - (level + 1) * decomp_base_log.0;
            let scaled_pt = PlaintextList::from_container(
                (0..input_sk_poly.polynomial_size().0)
                    .map(|i| *neg_sk_poly.get(i).0 << log_scale)
                    .collect::<Vec<u64>>(),
            );
            encrypt_glwe_ciphertext(
                output_glwe_sk,
                &mut glwe,
                &scaled_pt,
                Gaussian::from_dispersion_parameter(noise_parameters, 0.0),
                generator,
            );
        }
    }

    glwe_keyswitch_key
}

fn standard_keyswitch_glwe_ciphertext<InputCont, OutputCont>(
    glwe_keyswitch_key: &GlweKeyswitchKey,
    input_glwe_ciphertext: &GlweCiphertext<InputCont>,
    output_glwe_ciphertext: &mut GlweCiphertext<OutputCont>,
) where
    InputCont: Container<Element = u64>,
    OutputCont: ContainerMut<Element = u64>,
{
    let polynomial_size = glwe_keyswitch_key.polynomial_size();
    let output_glwe_size = glwe_keyswitch_key.output_glwe_dimension().to_glwe_size();
    let decomp_base_log = glwe_keyswitch_key.decomp_base_log();
    let decomp_level = glwe_keyswitch_key.decomp_level_count();
    let ciphertext_modulus = glwe_keyswitch_key.ciphertext_modulus();

    output_glwe_ciphertext.as_mut().fill(0u64);
    let (input_mask, input_body) = input_glwe_ciphertext.get_mask_and_body();
    output_glwe_ciphertext
        .get_mut_body()
        .as_mut()
        .clone_from_slice(input_body.as_ref());

    let decomposer = SignedDecomposer::new(decomp_base_log, decomp_level);

    for (index, input_mask_poly) in input_mask.as_polynomial_list().iter().enumerate() {
        let mut input_decomp_poly_list =
            PolynomialList::new(0u64, polynomial_size, PolynomialCount(decomp_level.0));

        for (i, value) in input_mask_poly.iter().enumerate() {
            let decomposition_iter = decomposer.decompose(*value);
            for (j, decomp_val) in decomposition_iter.into_iter().enumerate() {
                *input_decomp_poly_list
                    .get_mut(j)
                    .as_mut()
                    .get_mut(i)
                    .unwrap() = decomp_val.value();
            }
        }

        let glev = glwe_keyswitch_key.glev_list(index);
        for (decomp_poly, glwe) in input_decomp_poly_list.iter().zip(glev.iter().rev()) {
            let mut buf =
                GlweCiphertext::new(0u64, output_glwe_size, polynomial_size, ciphertext_modulus);
            for (mut buf_poly, glwe_poly) in buf
                .as_mut_polynomial_list()
                .iter_mut()
                .zip(glwe.as_polynomial_list().iter())
            {
                polynomial_wrapping_mul(&mut buf_poly, &decomp_poly, &glwe_poly);
            }
            glwe_ciphertext_add_assign(output_glwe_ciphertext, &buf);
        }
    }
}

#[derive(Clone)]
pub struct AutomorphKey {
    ksk: GlweKeyswitchKey,
    glwe_dimension: GlweDimension,
    polynomial_size: PolynomialSize,
    auto_k: usize,
}

impl AutomorphKey {
    fn allocate(
        decomp_base_log: DecompositionBaseLog,
        decomp_level_count: DecompositionLevelCount,
        glwe_dimension: GlweDimension,
        polynomial_size: PolynomialSize,
        auto_k: usize,
    ) -> Self {
        Self {
            ksk: GlweKeyswitchKey::new(
                glwe_dimension,
                glwe_dimension,
                polynomial_size,
                decomp_base_log,
                decomp_level_count,
                CiphertextModulus::<u64>::new_native(),
            ),
            glwe_dimension,
            polynomial_size,
            auto_k,
        }
    }

    fn fill_with_automorph_key<G: ByteRandomGenerator>(
        &mut self,
        before_key: &mut GlweSecretKey<Vec<u64>>,
        after_key: &GlweSecretKey<Vec<u64>>,
        k: usize,
        noise_parameters: impl DispersionParameter,
        generator: &mut EncryptionRandomGenerator<G>,
    ) {
        let mut before_poly_list = PolynomialList::new(
            0u64,
            self.polynomial_size,
            PolynomialCount(self.glwe_dimension.0),
        );
        for (mut before_poly, after_poly) in before_poly_list
            .iter_mut()
            .zip(after_key.as_polynomial_list().iter())
        {
            let out = eval_x_k(after_poly.as_view(), k);
            before_poly.as_mut().clone_from_slice(out.as_ref());
        }
        *before_key =
            GlweSecretKey::from_container(before_poly_list.into_container(), self.polynomial_size);

        self.ksk = allocate_and_generate_new_glwe_keyswitch_key(
            before_key,
            after_key,
            self.ksk.decomp_base_log(),
            self.ksk.decomp_level_count(),
            noise_parameters,
            CiphertextModulus::<u64>::new_native(),
            generator,
        );
        self.auto_k = k;
    }

    pub(crate) fn auto<InputCont, OutputCont>(
        &self,
        after: &mut GlweCiphertext<OutputCont>,
        before: &GlweCiphertext<InputCont>,
    ) where
        InputCont: Container<Element = u64>,
        OutputCont: ContainerMut<Element = u64>,
    {
        let mut before_power = GlweCiphertextOwned::new(
            0u64,
            before.glwe_size(),
            before.polynomial_size(),
            before.ciphertext_modulus(),
        );
        for (mut poly_power, poly) in before_power
            .as_mut_polynomial_list()
            .iter_mut()
            .zip(before.as_polynomial_list().iter())
        {
            poly_power
                .as_mut()
                .clone_from_slice(eval_x_k(poly, self.auto_k).as_ref());
        }

        standard_keyswitch_glwe_ciphertext(&self.ksk, &before_power, after);
    }
}

pub fn gen_all_auto_keys<G: ByteRandomGenerator>(
    decomp_base_log: DecompositionBaseLog,
    decomp_level: DecompositionLevelCount,
    glwe_secret_key: &GlweSecretKey<Vec<u64>>,
    noise_parameters: impl DispersionParameter,
    generator: &mut EncryptionRandomGenerator<G>,
) -> HashMap<usize, AutomorphKey> {
    let glwe_dimension = glwe_secret_key.glwe_dimension();
    let polynomial_size = glwe_secret_key.polynomial_size();

    let mut keys = HashMap::new();
    for i in 1..=(polynomial_size.0).ilog2() as usize {
        let k = polynomial_size.0 / (1 << (i - 1)) + 1;
        let mut auto_key = AutomorphKey::allocate(
            decomp_base_log,
            decomp_level,
            glwe_dimension,
            polynomial_size,
            i,
        );
        let mut before_key = glwe_secret_key.clone();
        auto_key.fill_with_automorph_key(
            &mut before_key,
            glwe_secret_key,
            k,
            noise_parameters,
            generator,
        );
        keys.insert(k, auto_key);
    }

    keys
}
