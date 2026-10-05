use crate::arithmetic::{self, ArithmeticCtx};
use crate::cbs;
use crate::lut_eval;
use crate::pbs;
use crate::utils;
use rayon::join;
use rayon::prelude::*;
use tfhe::core_crypto::entities::GlweCiphertext;
use tfhe::core_crypto::prelude::*;

pub fn num_blocks(m: usize, polynomial_size: usize) -> usize {
    m.div_ceil(polynomial_size / 2)
}

#[allow(clippy::too_many_arguments)]
pub fn extract_subsequence(
    ct_text_parts: &[GlweCiphertextList<Vec<u64>>],
    ct_text_shifted_parts: &[GlweCiphertextList<Vec<u64>>],
    i_parts: &[LweCiphertext<Vec<u64>>],
    m: usize,
    n_parts: usize,
    top_part_bits: usize,
    pfpksk: &LwePrivateFunctionalPackingKeyswitchKeyList<Vec<u64>>,
    cbs_base_log: DecompositionBaseLog,
    cbs_level: DecompositionLevelCount,
    ctx: &ArithmeticCtx<'_>,
) -> Vec<LweCiphertext<Vec<u64>>> {
    let parts_per_char = ct_text_parts.len();
    assert_eq!(ct_text_shifted_parts.len(), parts_per_char);
    assert!(parts_per_char > 0);

    let half_n = ctx.polynomial_size.0 / 2;
    let b = num_blocks(m, ctx.polynomial_size.0);

    // Offsets k*half_n for k = 1..=b, i.e. the b additional index values
    // (besides i_parts itself) shared between the A_j (offset j) and B_j
    // (offset j+1) series. k*half_n is a known public constant, so each is
    // computed via `add_plain_offset_lwe_vec` (no ciphertext encryption of
    // the offset, no RNG), independently and in parallel.
    let mut offset_index_parts: Vec<Vec<LweCiphertext<Vec<u64>>>> = Vec::with_capacity(b + 1);
    offset_index_parts.push(i_parts.to_vec());
    let offsets: Vec<Vec<LweCiphertext<Vec<u64>>>> = (1..=b)
        .into_par_iter()
        .map(|k| {
            let summed = arithmetic::add_plain_offset_lwe_vec(i_parts, (k * half_n) as u64, ctx);
            summed[..n_parts].to_vec()
        })
        .collect();
    offset_index_parts.extend(offsets);

    // Step 1: GGSW generation for each of the b+1 offsets, parallel across
    // offsets and, within each offset, across the n_parts radix nibbles.
    let ggsw_lists: Vec<GgswCiphertextList<Vec<u64>>> = offset_index_parts
        .par_iter()
        .map(|parts| {
            let ggsw_parts: Vec<_> = parts
                .par_iter()
                .enumerate()
                .map(|(idx, ct)| {
                    let bits = if idx == n_parts - 1 { top_part_bits } else { 2 };
                    cbs::cbs_int_4to2(ct, ctx.ksk, ctx.bsk, pfpksk, bits)
                })
                .collect();
            utils::concat_ggsw_lists(&ggsw_parts)
        })
        .collect();

    // Step 2: CMux-tree evaluation for one part-series. ggsw_lists[j] feeds
    // A_j via ct_text, ggsw_lists[j+1] feeds B_j via ct_text_shifted; the two
    // series run in parallel (rayon::join), each itself parallel across its b
    // blocks. This is repeated once per character-part series (parts_per_char
    // of them); the index-dependent GGSW generation above (Step 1) is shared.
    let compute_blocks = |ct_text: &GlweCiphertextList<Vec<u64>>,
                          ct_text_shifted: &GlweCiphertextList<Vec<u64>>|
     -> (Vec<GlweCiphertext<Vec<u64>>>, Vec<GlweCiphertext<Vec<u64>>>) {
        join(
            || {
                (0..b)
                    .into_par_iter()
                    .map(|j| {
                        lut_eval::private_lut_eval_n_to_4_rotated_rlwe(
                            ct_text,
                            &ggsw_lists[j],
                            ctx.polynomial_size,
                            cbs_base_log,
                            cbs_level,
                            ctx.ciphertext_modulus,
                        )
                    })
                    .collect()
            },
            || {
                (0..b)
                    .into_par_iter()
                    .map(|j| {
                        lut_eval::private_lut_eval_n_to_4_rotated_rlwe(
                            ct_text_shifted,
                            &ggsw_lists[j + 1],
                            ctx.polynomial_size,
                            cbs_base_log,
                            cbs_level,
                            ctx.ciphertext_modulus,
                        )
                    })
                    .collect()
            },
        )
    };

    let (a_blocks_0, b_blocks_0) = compute_blocks(&ct_text_parts[0], &ct_text_shifted_parts[0]);

    // Nega check: computed once from part 0's A_0 last coefficient, exactly
    // as in the b=1 case. This reflects the negacyclic wrap-around of the
    // encrypted index/rotation, not the text content, so it is shared across
    // every character-part series below.
    let mut ct_last_char = LweCiphertext::new(
        0u64,
        ctx.lwe_secret_key.lwe_dimension().to_lwe_size(),
        ctx.ciphertext_modulus,
    );
    extract_lwe_sample_from_glwe_ciphertext(
        &a_blocks_0[0],
        &mut ct_last_char,
        MonomialDegree(half_n - 1),
    );
    let mut selector_pm1 = LweCiphertext::new(0u64, ct_last_char.lwe_size(), ctx.ciphertext_modulus);
    pbs::pbs_selector_pm1_from_last_char(
        &ct_last_char,
        ctx.ksk,
        ctx.bsk,
        ctx.lwe_secret_key_after_ks,
        ctx.polynomial_size,
        ctx.glwe_size,
        ctx.ciphertext_modulus,
        ctx.delta,
        &mut selector_pm1,
    );
    let sel = cbs::cbs_int_4to1(&selector_pm1, ctx.ksk, ctx.bsk, pfpksk, 3);
    let fourier_sel = {
        let mut fourier = FourierGgswCiphertext::new(
            ctx.glwe_size,
            ctx.polynomial_size,
            cbs_base_log,
            cbs_level,
        );
        convert_standard_ggsw_ciphertext_to_fourier(&sel, &mut fourier);
        fourier
    };

    // Step 3: for each block j, select A_j/B_j according to r_j = s xor (j mod 2)
    // -- implemented by swapping the cmux argument order instead of negating s
    // -- then extract its first half_n LWE samples, truncated to m. Parallel
    // across blocks.
    let finalize = |a_blocks: Vec<GlweCiphertext<Vec<u64>>>,
                    b_blocks: Vec<GlweCiphertext<Vec<u64>>>|
     -> Vec<LweCiphertext<Vec<u64>>> {
        let blocks: Vec<Vec<LweCiphertext<Vec<u64>>>> = a_blocks
            .into_iter()
            .zip(b_blocks)
            .enumerate()
            .collect::<Vec<_>>()
            .into_par_iter()
            .map(|(j, (mut a_j, mut b_j))| {
                let c_j = if j % 2 == 0 {
                    cmux_assign(&mut a_j, &mut b_j, &fourier_sel);
                    a_j
                } else {
                    cmux_assign(&mut b_j, &mut a_j, &fourier_sel);
                    b_j
                };
                (0..half_n)
                    .map(|pos| {
                        let mut out = LweCiphertext::new(
                            0u64,
                            ctx.lwe_secret_key.lwe_dimension().to_lwe_size(),
                            ctx.ciphertext_modulus,
                        );
                        extract_lwe_sample_from_glwe_ciphertext(&c_j, &mut out, MonomialDegree(pos));
                        out
                    })
                    .collect()
            })
            .collect();

        blocks.into_iter().flatten().take(m).collect()
    };

    // Per-character-part results (part 0 is already computed above).
    let mut per_part: Vec<Vec<LweCiphertext<Vec<u64>>>> = Vec::with_capacity(parts_per_char);
    per_part.push(finalize(a_blocks_0, b_blocks_0));
    if parts_per_char > 1 {
        let rest: Vec<Vec<LweCiphertext<Vec<u64>>>> = (1..parts_per_char)
            .into_par_iter()
            .map(|k| {
                let (a_k, b_k) = compute_blocks(&ct_text_parts[k], &ct_text_shifted_parts[k]);
                finalize(a_k, b_k)
            })
            .collect();
        per_part.extend(rest);
    }

    // Interleave char-major, part-minor (MSB -> LSB), matching the pattern's
    // ciphertext ordering expected by `sgn`.
    let mut interleaved = Vec::with_capacity(m * parts_per_char);
    for char_idx in 0..m {
        for part_k in &per_part {
            interleaved.push(part_k[char_idx].clone());
        }
    }
    interleaved
}
