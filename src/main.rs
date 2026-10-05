use clap::Parser;
use rand::{Rng, SeedableRng};
use std::path::PathBuf;
use tfhe::core_crypto::algorithms::{
    allocate_and_generate_new_circuit_bootstrap_lwe_pfpksk_list, encrypt_lwe_ciphertext,
};
use tfhe::core_crypto::prelude::*;
use tfhe_pattern::{
    automorphism::gen_all_auto_keys,
    fhe_find,
    keygen::keygen_pbs,
    param::{self, FheKeys, FheParams},
};

#[derive(Clone, Copy, Debug, clap::ValueEnum)]
enum Alphabet {
    Dna,
    Ascii,
}

impl Alphabet {
    fn parts_per_char(self) -> usize {
        match self {
            Alphabet::Dna => 1,
            Alphabet::Ascii => 4,
        }
    }
}

/// Decompose a raw symbol code into `parts_per_char` radix-4 parts, ordered
/// MSB -> LSB, each encoded as (2-bit value + 1) so that normal symbols stay
/// in 1..=4 (0 is reserved for pattern padding, 5 for the text sentinel `$`).
/// Smallest `log_n` such that `2^log_n >= x` (for `x >= 1`).
fn ceil_log2(x: u64) -> usize {
    assert!(x >= 1, "❗️--n must be >= 1");
    if x == 1 {
        0
    } else {
        (64 - (x - 1).leading_zeros()) as usize
    }
}

fn encode_symbol(raw: u64, parts_per_char: usize) -> Vec<u64> {
    if parts_per_char == 1 {
        vec![raw]
    } else {
        (0..parts_per_char)
            .map(|k| {
                let shift = 2 * (parts_per_char - 1 - k);
                ((raw >> shift) & 0x3) + 1
            })
            .collect()
    }
}

/// Load raw bytes from an inline string (as-is, UTF-8 bytes) or a file
/// (as-is, including any trailing newline), or `None` if neither was given.
/// `--text`/`--text-file` (and `--pattern`/`--pattern-file`) are mutually
/// exclusive via clap's `conflicts_with`, so at most one of the two is set.
fn load_explicit_bytes(inline: &Option<String>, file: &Option<PathBuf>) -> Option<Vec<u8>> {
    if let Some(s) = inline {
        Some(s.clone().into_bytes())
    } else {
        file.as_ref().map(|path| {
            std::fs::read(path)
                .unwrap_or_else(|e| panic!("❗️failed to read {}: {}", path.display(), e))
        })
    }
}

/// Decode raw bytes into the per-character raw codes used elsewhere in this
/// program (DNA: 1..=4; ASCII: 0..=127), validating that every byte belongs
/// to the selected alphabet.
fn decode_explicit_bytes(bytes: &[u8], alphabet: Alphabet) -> Vec<u64> {
    match alphabet {
        Alphabet::Ascii => bytes
            .iter()
            .map(|&b| {
                if b > 0x7f {
                    panic!(
                        "❗️--alphabet ascii requires 7-bit ASCII input (0x00..=0x7f), found byte 0x{:02x}",
                        b
                    );
                }
                b as u64
            })
            .collect(),
        Alphabet::Dna => bytes
            .iter()
            .map(|&b| match b.to_ascii_uppercase() {
                b'A' => 1,
                b'C' => 2,
                b'G' => 3,
                b'T' => 4,
                other => panic!(
                    "❗️--alphabet dna requires only A/C/G/T characters, found byte 0x{:02x}",
                    other
                ),
            })
            .collect(),
    }
}

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Character alphabet: "dna" (A/C/G/T, parts_per_char=1) or "ascii"
    /// (7-bit ASCII, parts_per_char=4)
    #[arg(long, value_enum, default_value = "ascii")]
    alphabet: Alphabet,

    /// RNG seed (random if None)
    #[arg(long)]
    seed: Option<u64>,

    /// length of text T after padding, as log2(n); n = 2^log_n (log_n >= 11).
    /// Mutually exclusive with --n. If neither is given, the smallest
    /// capacity satisfying the NFBE safety margin is chosen automatically.
    #[arg(long)]
    log_n: Option<usize>,

    /// minimum padded text capacity; rounded up to the next power of two.
    /// Mutually exclusive with --log-n.
    #[arg(long, conflicts_with = "log_n")]
    n: Option<u64>,

    /// length of text T before padding (n0 <= n); ignored if --text/--text-file is given
    #[arg(long, default_value_t = 1000)]
    n0: usize,

    /// length of pattern P after padding; defaults to m0 (no wildcard
    /// padding) if not given
    #[arg(long)]
    m: Option<usize>,

    /// length of pattern P before padding (m0 <= m); ignored if --pattern/--pattern-file is given
    #[arg(long, default_value_t = 20)]
    m0: usize,

    /// Explicit text content (n0 is derived from its length); conflicts with --text-file
    #[arg(long, conflicts_with = "text_file")]
    text: Option<String>,

    /// Read the text content from a file, as-is (n0 is derived from its length)
    #[arg(long)]
    text_file: Option<PathBuf>,

    /// Explicit pattern content (m0 is derived from its length); conflicts with --pattern-file
    #[arg(long, conflicts_with = "pattern_file")]
    pattern: Option<String>,

    /// Read the pattern content from a file, as-is (m0 is derived from its length)
    #[arg(long)]
    pattern_file: Option<PathBuf>,

    /// T[solution: solution + m0] == P (random if None, random pattern if -1)
    #[arg(long)]
    solution: Option<i64>,

    /// Measure only suffix-array construction
    #[arg(long, default_value_t = false)]
    sa_only: bool,
}

#[allow(non_snake_case)]
fn main() {
    let args = Args::parse();

    /* Client Side */
    let mut rng = rand::rng();
    let seed = match args.seed {
        Some(s) => s,
        None => rng.random::<u64>(),
    };
    let mut rng = rand::rngs::StdRng::seed_from_u64(seed);
    // Explicit text/pattern (--text/--text-file, --pattern/--pattern-file)
    // take priority over the random-generation mode; when given, n0/m0 are
    // derived from their length instead of --n0/--m0.
    let text_bytes = load_explicit_bytes(&args.text, &args.text_file);
    let pattern_bytes = load_explicit_bytes(&args.pattern, &args.pattern_file);
    let n0: usize = text_bytes.as_ref().map_or(args.n0, |b| b.len());
    let m0: usize = pattern_bytes.as_ref().map_or(args.m0, |b| b.len());
    // Default to no wildcard padding (m = m0) unless the user explicitly
    // asks for extra padded positions via --m.
    let m: usize = args.m.unwrap_or(m0);
    let solution: Option<i64> = args.solution;
    let lwe_dimension = param::lwe_dimension();
    let glwe_dimension = param::glwe_dimension();
    let glwe_size = glwe_dimension.to_glwe_size();
    let polynomial_size = param::polynomial_size();
    let N: usize = polynomial_size.0;
    let lwe_modular_std_dev = param::lwe_modular_std_dev();
    let glwe_modular_std_dev = param::glwe_modular_std_dev();
    let pbs_base_log = param::pbs_base_log();
    let pbs_level = param::pbs_level();
    let ks_base_log = param::ks_base_log();
    let ks_level = param::ks_level();
    let cbs_base_log = param::cbs_base_log();
    let cbs_level = param::cbs_level();
    let auto_base_log = param::auto_base_log();
    let auto_level = param::auto_level();
    let ss_base_log = param::ss_base_log();
    let ss_level = param::ss_level();
    let ciphertext_modulus = param::ciphertext_modulus();
    // Independent of log_n/n, so it can be computed before log_n is resolved
    // (and reused below instead of recomputing it from scratch).
    let n_blocks = tfhe_pattern::nfbe::num_blocks(m, N);
    // log_n/n resolution: --log-n fixes log_n exactly; --n rounds up to the
    // next power of two (log_n = ceil(log2(n))); with neither, the smallest
    // capacity satisfying the NFBE safety margin (n > n0 + n_blocks*(N/2))
    // is chosen automatically. All three are clamped to the log_n >= 11
    // floor enforced just below.
    let log_n: usize = match (args.log_n, args.n) {
        (Some(log_n_arg), None) => log_n_arg,
        (None, Some(n_arg)) => ceil_log2(n_arg).max(11),
        (None, None) => {
            let min_capacity = n0 as u64 + (n_blocks * (N / 2)) as u64 + 1;
            ceil_log2(min_capacity).max(11)
        }
        (Some(_), Some(_)) => unreachable!("clap enforces --log-n/--n are mutually exclusive"),
    };
    let n: u64 = 2u64.pow(log_n as u32);
    const INDEX_PART_BITS: usize = 2;
    let n_parts = log_n.div_ceil(INDEX_PART_BITS);
    let top_part_bits = {
        let rem = log_n % INDEX_PART_BITS;
        if rem == 0 { INDEX_PART_BITS } else { rem }
    };
    let mut boxed_seeder = new_seeder();
    let seeder = boxed_seeder.as_mut();
    let mut secret_generator = SecretRandomGenerator::<DefaultRandomGenerator>::new(seeder.seed());
    let mut encryption_generator =
        EncryptionRandomGenerator::<DefaultRandomGenerator>::new(seeder.seed(), seeder);
    let delta = 1u64 << 60;
    if log_n < 11 {
        panic!("❗️log_n >= 11 is required");
    }
    if n < n0 as u64 || m < m0 || m0 > n0 {
        panic!("❗️n0 <= n, m0 <= m, and m0 <= n0 are required");
    }
    if (n0 as u64 + (n_blocks * (N / 2)) as u64) >= n {
        panic!(
            "❗️n0 + b*(N/2) < n is required (N = {}, b = {})",
            N, n_blocks
        );
    }
    // `solution` only drives the random self-consistent pattern generation
    // below; it is meaningless (and ignored) once an explicit pattern is
    // given, since that pattern need not be a substring copied from T.
    if pattern_bytes.is_none()
        && let Some(s) = solution
        && s >= 0
        && (s as u64 + m0 as u64 > n0 as u64)
    {
        panic!("❗️solution + m0 <= n0 is required");
    }
    let solution_val: u64 = if pattern_bytes.is_some() {
        0 // unused: pattern is explicit, not copied from T at an offset
    } else {
        match solution {
            Some(s) if s >= 0 => s as u64,
            Some(-1) => 0, // dummy, not used
            None => rng.random_range(0..=(n0 as u64 - m0 as u64)),
            _ => panic!("solution must be >= -1"),
        }
    };
    println!("log_n = {}", log_n);
    println!("n     = {}", n);
    println!("n0    = {}", n0);
    println!("m     = {}", m);
    println!("m0    = {}", m0);
    println!("N     = {}", N);
    println!(
        "sol   = {}",
        if pattern_bytes.is_some() {
            "N/A (explicit pattern)".to_string()
        } else {
            match solution {
                Some(s) => s.to_string(),
                None => solution_val.to_string(),
            }
        }
    );
    println!("Seed  = {}", seed);
    println!("n_parts       = {}", n_parts);
    println!("top_part_bits = {}", top_part_bits);
    let parts_per_char = args.alphabet.parts_per_char();
    println!("alphabet      = {:?}", args.alphabet);
    println!("parts_per_char = {}", parts_per_char);
    // Raw per-character codes: DNA uses 1..=4 (A/C/G/T) directly as the
    // "parts_per_char=1" case of encode_symbol; ASCII uses 0..=127. Padding
    // (pattern, value 0) and the text sentinel `$` (value 5) are injected
    // uniformly across all parts_per_char parts below, not decomposed via
    // encode_symbol.
    let t_raw: Vec<u64> = match &text_bytes {
        Some(bytes) => decode_explicit_bytes(bytes, args.alphabet),
        None => (0..n0)
            .map(|_| match args.alphabet {
                Alphabet::Dna => rng.random_range(1..=4),
                Alphabet::Ascii => rng.random_range(0..128),
            })
            .collect(),
    };
    let pattern_raw: Vec<u64> = match &pattern_bytes {
        Some(bytes) => decode_explicit_bytes(bytes, args.alphabet),
        None if solution == Some(-1) => (0..m0)
            .map(|_| match args.alphabet {
                Alphabet::Dna => rng.random_range(1..=4),
                Alphabet::Ascii => rng.random_range(0..128),
            })
            .collect(),
        None => (0..m0).map(|i| t_raw[solution_val as usize + i]).collect(),
    };
    // parts[k][i] = part k (MSB -> LSB) of the character at position i;
    // positions beyond the real content are filled with `sentinel`.
    let build_char_parts = |raw: &[u64], total_len: usize, sentinel: u64| -> Vec<Vec<u64>> {
        let mut parts = vec![vec![sentinel; total_len]; parts_per_char];
        for (i, &v) in raw.iter().enumerate() {
            let sym = encode_symbol(v, parts_per_char);
            for (k, part) in parts.iter_mut().enumerate() {
                part[i] = sym[k];
            }
        }
        parts
    };
    let t_parts = build_char_parts(&t_raw, n as usize, 5);
    let pattern_parts = build_char_parts(&pattern_raw, m, 0);
    let text_preview = |raw: &[u64]| -> String {
        match args.alphabet {
            Alphabet::Dna => raw
                .iter()
                .map(|&v| match v {
                    1 => "A",
                    2 => "C",
                    3 => "G",
                    4 => "T",
                    _ => "?",
                })
                .collect::<Vec<_>>()
                .join(""),
            Alphabet::Ascii => raw.iter().map(|&v| (v as u8 as char).to_string()).collect(),
        }
    };
    println!(
        "T         : {}",
        text_preview(&t_raw[..t_raw.len().min(30)])
    );
    println!(
        "P         : {}",
        text_preview(&pattern_raw[..pattern_raw.len().min(30)])
    );
    println!("========== Suffix Array Construction ==========");
    let start_sa = std::time::Instant::now();
    // Suffix array of the real text T[0..n0-1] only: out-of-range comparisons
    // use usize::MAX as sentinel, which (like the dummy symbol 5 used to fill
    // T[n0..n-1]) always compares greater than any real symbol (1..=4), so the
    // relative order of real suffixes is identical to sorting the full padded
    // array of length n.
    let sa = {
        let n = n0;
        let mut sa: Vec<usize> = (0..n).collect();
        let mut rank: Vec<usize> = t_raw.iter().map(|&v| v as usize).collect();
        let mut tmp = vec![0usize; n];
        let mut k = 1usize;
        while {
            sa.sort_by_key(|&i| {
                let r1 = rank[i];
                let r2 = if i + k < n { rank[i + k] } else { usize::MAX };
                (r1, r2)
            });
            tmp[sa[0]] = 0;
            for i in 1..n {
                let a = sa[i - 1];
                let b = sa[i];
                let left = (rank[a], if a + k < n { rank[a + k] } else { usize::MAX });
                let right = (rank[b], if b + k < n { rank[b + k] } else { usize::MAX });
                tmp[b] = tmp[a] + if left != right { 1 } else { 0 };
            }
            rank.copy_from_slice(&tmp);
            k <<= 1;
            rank[sa[n - 1]] != n - 1
        } {}
        sa
    };
    let elapsed_sa = start_sa.elapsed();
    println!(
        "🧮 {:.3} [s] in suffix array construction",
        elapsed_sa.as_secs_f64()
    );
    if args.sa_only {
        println!("SA length = {}", sa.len());
        return;
    }
    // Pad the SA out to length n for encryption: every rank beyond n0-1
    // corresponds to a dummy suffix, and is mapped to the fixed, safe index
    // n0 (the start of the dummy region) rather than its true (and otherwise
    // unbounded, up to n-1) text position. This guarantees that reading
    // T[i:i+b*(N/2)-1] for any dummy rank stays within n0+b*(N/2)<=n and
    // always yields the all-dummy sequence, instead of risking a wrap-around
    // read of real text near the start of T.
    let sa_full: Vec<usize> = (0..n as usize)
        .map(|rank| if rank < n0 { sa[rank] } else { n0 })
        .collect();
    let lut_parts = {
        let len = sa_full.len();
        let mut parts = vec![vec![0u64; len]; n_parts];
        for i in 0..len {
            let v = sa_full[i] as u64;
            for part in 0..n_parts {
                parts[part][i] = ((v >> (part * 2)) & 0x3) as u64;
            }
        }
        parts
    };
    // Same shift-by-N/2 construction as T, applied independently to each
    // character part.
    let t_shifted_parts: Vec<Vec<u64>> = t_parts
        .iter()
        .map(|part| {
            let mut shifted = vec![5u64; n as usize];
            for i in N / 2..n as usize {
                let idx = i - N / 2;
                if idx < n as usize {
                    shifted[i] = part[idx];
                }
            }
            shifted
        })
        .collect();
    println!("========== Key Generation ==========");
    let start_keygen = std::time::Instant::now();
    let (lwe_secret_key, glwe_secret_key, lwe_secret_key_after_ks, bsk, ksk) = keygen_pbs(
        lwe_dimension,
        glwe_dimension,
        polynomial_size,
        lwe_modular_std_dev,
        glwe_modular_std_dev,
        pbs_base_log,
        pbs_level,
        ks_base_log,
        ks_level,
        &mut secret_generator,
        &mut encryption_generator,
    );
    let auto_keys = gen_all_auto_keys(
        auto_base_log,
        auto_level,
        &glwe_secret_key,
        glwe_modular_std_dev,
        &mut encryption_generator,
    );
    let pfpksk = allocate_and_generate_new_circuit_bootstrap_lwe_pfpksk_list(
        &lwe_secret_key,
        &glwe_secret_key,
        ss_base_log,
        ss_level,
        Gaussian::from_dispersion_parameter(glwe_modular_std_dev, 0.0),
        ciphertext_modulus,
        &mut encryption_generator,
    );
    let mut ct_text: Vec<GlweCiphertextList<Vec<u64>>> = (0..parts_per_char)
        .map(|_| {
            GlweCiphertextList::new(
                0u64,
                glwe_size,
                polynomial_size,
                GlweCiphertextCount((n / N as u64) as usize),
                ciphertext_modulus,
            )
        })
        .collect();
    let mut ct_text_shifted: Vec<GlweCiphertextList<Vec<u64>>> = (0..parts_per_char)
        .map(|_| {
            GlweCiphertextList::new(
                0u64,
                glwe_size,
                polynomial_size,
                GlweCiphertextCount((n / N as u64) as usize),
                ciphertext_modulus,
            )
        })
        .collect();
    let mut ct_sa: Vec<GlweCiphertextList<Vec<u64>>> = (0..n_parts)
        .map(|_| {
            GlweCiphertextList::new(
                0u64,
                glwe_size,
                polynomial_size,
                GlweCiphertextCount((n / N as u64) as usize),
                ciphertext_modulus,
            )
        })
        .collect();
    let noise_distribution = Gaussian::from_dispersion_parameter(glwe_modular_std_dev, 0.0);
    for part in 0..parts_per_char {
        for (glwe_idx, mut glwe) in ct_text[part].iter_mut().enumerate() {
            let mut plaintext_poly = vec![0u64; polynomial_size.0];
            for (i, val) in plaintext_poly.iter_mut().enumerate() {
                *val = t_parts[part][glwe_idx * N + i] * delta;
            }
            encrypt_glwe_ciphertext(
                &glwe_secret_key,
                &mut glwe,
                &PlaintextList::from_container(plaintext_poly),
                noise_distribution,
                &mut encryption_generator,
            );
        }
        for (glwe_idx, mut glwe) in ct_text_shifted[part].iter_mut().enumerate() {
            let mut plaintext_poly = vec![0u64; polynomial_size.0];
            for (i, val) in plaintext_poly.iter_mut().enumerate() {
                *val = t_shifted_parts[part][glwe_idx * N + i] * delta;
            }
            encrypt_glwe_ciphertext(
                &glwe_secret_key,
                &mut glwe,
                &PlaintextList::from_container(plaintext_poly),
                noise_distribution,
                &mut encryption_generator,
            );
        }
    }
    for part in 0..n_parts {
        for (glwe_idx, mut cur_glwe) in ct_sa[part].iter_mut().enumerate() {
            let mut plaintext_poly = vec![0u64; polynomial_size.0];
            for (i, val) in plaintext_poly.iter_mut().enumerate() {
                *val = (lut_parts[part][glwe_idx * N + i] as u64) * delta;
            }
            encrypt_glwe_ciphertext(
                &glwe_secret_key,
                &mut cur_glwe,
                &PlaintextList::from_container(plaintext_poly),
                noise_distribution,
                &mut encryption_generator,
            );
        }
    }
    // Flatten pattern_parts (parts_per_char x m) into char-major, part-minor
    // (MSB -> LSB) order: char0_part0, char0_part1, ..., char1_part0, ...
    // matching the ordering `nfbe::extract_subsequence` produces for T.
    let mut pattern_flat: Vec<u64> = Vec::with_capacity(m * parts_per_char);
    for char_idx in 0..m {
        for part in &pattern_parts {
            pattern_flat.push(part[char_idx]);
        }
    }
    let ct_pattern: Vec<LweCiphertext<Vec<u64>>> = pattern_flat
        .iter()
        .map(|&v| {
            let mut ct = LweCiphertext::new(
                0u64,
                lwe_secret_key.lwe_dimension().to_lwe_size(),
                ciphertext_modulus,
            );
            encrypt_lwe_ciphertext(
                &lwe_secret_key,
                &mut ct,
                Plaintext(v * delta),
                noise_distribution,
                &mut encryption_generator,
            );
            ct
        })
        .collect();
    let keys = FheKeys {
        lwe_secret_key: lwe_secret_key.clone(),
        glwe_secret_key,
        lwe_secret_key_after_ks,
        bsk,
        ksk,
        pfpksk,
        auto_keys,
    };
    let fhe_params = FheParams {
        log_n,
        n,
        m,
        n_parts,
        top_part_bits,
        polynomial_size,
        glwe_size,
        glwe_modular_std_dev,
        cbs_base_log,
        cbs_level,
        ciphertext_modulus,
        encryption_generator,
        delta,
    };
    let elapsed_keygen = start_keygen.elapsed();
    println!(
        "🔑 {:.3} [s] in key generation",
        elapsed_keygen.as_secs_f64()
    );

    /* Server Side */
    let start_find = std::time::Instant::now();
    println!("========== FHE find ==========");
    let result = fhe_find(
        ct_text,
        ct_text_shifted,
        ct_sa,
        ct_pattern,
        keys,
        fhe_params,
    );
    let elapsed_find = start_find.elapsed();
    println!("🔍 {:.3} [s] in FHE find", elapsed_find.as_secs_f64());

    /* Client Side (we omit the part of decryption request of sol + r to the Owner, where r is a random value)  */
    let mut sol_recon: u64 = 0;
    for (i, ct) in result.index.iter().enumerate() {
        let decrypted = decrypt_lwe_ciphertext(&lwe_secret_key, ct).0;
        let part = ((decrypted + (delta / 2)) / delta) as u64 & 0x3;
        sol_recon |= part << (2 * i);
    }
    let found_dec = decrypt_lwe_ciphertext(&lwe_secret_key, &result.found).0;
    let found = ((found_dec + (delta / 2)) / delta) as u64 & 0x1;

    let sol_pos = sol_recon as usize;
    if found == 1 {
        if sol_pos + m0 <= t_raw.len() {
            let t_substring: Vec<u64> = t_raw[sol_pos..sol_pos + m0].to_vec();
            let match_found = t_substring == pattern_raw[..m0].to_vec();
            if match_found {
                println!("✅ Match: T[{}:{}] == P", sol_pos, sol_pos + m0);
            } else {
                println!(
                    "❌ No Match (found=1 but decoded index mismatch: {})",
                    sol_recon
                );
            }
        } else {
            println!(
                "❌ No Match (found=1 but decoded index out of range: {})",
                sol_recon
            );
        }
    } else {
        println!("❌ No Match");
    }
}
