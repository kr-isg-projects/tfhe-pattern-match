use tfhe::core_crypto::prelude::*;
use tfhe::core_crypto::algorithms::encrypt_lwe_ciphertext;
use tfhe_pattern::{
    fhe_find,
    param::FheKeys,
    param::FheParams,
    param::HYBRID_BASE_64,
};
use refined_tfhe_lhe::{
    keygen_pbs,
    gen_all_auto_keys,
    generate_scheme_switching_key,
};
use rand::{Rng, SeedableRng};
use clap::Parser;


#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// RNG seed (random if None)
    #[arg(long)]
    seed: Option<u64>,

    /// length of text T after padding (log-n >= 12)
    #[arg(long, default_value_t = 12)]
    log_n: usize,

    /// length of text T before padding (n0 <= n)
    #[arg(long, default_value_t = 1000)]
    n0: usize,

    /// length of pattern P after padding (m <= 2048)
    #[arg(long, default_value_t = 20)]
    m: usize,

    /// length of pattern P before padding (m0 <= m)
    #[arg(long, default_value_t = 20)]
    m0: usize,

    /// T[solution: solution + m0] == P (random if None, random pattern if -1)
    #[arg(long)]
    solution: Option<i64>,
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
    let log_n = args.log_n;
    let n: u64 = 2u64.pow(log_n as u32);
    let n0: usize = args.n0;
    let m: usize = args.m;
    let m0: usize = args.m0;
    let solution: Option<i64> = args.solution;
    let solution_val: u64 = match solution {
        Some(s) if s >= 0 => s as u64,
        Some(-1) => 0, // dummy, not used
        None => rng.random_range(0..=(n0 as u64 - m0 as u64)),
        _ => panic!("solution must be >= -1"),
    };
    let N: usize = 4096;
    let n_nibbles = (log_n as f64 / 4.0).ceil() as usize;
    let n_rems = if log_n % 4 == 0 { 4 } else { log_n % 4 };
    let param = *HYBRID_BASE_64;
    let lwe_dimension = param.lwe_dimension();
    let glwe_dimension = param.glwe_dimension();
    let glwe_size = glwe_dimension.to_glwe_size();
    let polynomial_size = param.polynomial_size();
    let lwe_modular_std_dev = param.lwe_modular_std_dev();
    let glwe_modular_std_dev = param.glwe_modular_std_dev();
    let pbs_base_log = param.pbs_base_log();
    let pbs_level = param.pbs_level();
    let ks_base_log = param.ks_base_log();
    let ks_level = param.ks_level();
    let cbs_base_log = param.cbs_base_log();
    let cbs_level = param.cbs_level();
    let auto_base_log = param.auto_base_log();
    let auto_level = param.auto_level();
    let auto_fft_type = param.fft_type_auto();
    let ss_base_log = param.ss_base_log();
    let ss_level = param.ss_level();
    let ciphertext_modulus = param.ciphertext_modulus();
    let mut boxed_seeder = new_seeder();
    let seeder = boxed_seeder.as_mut();
    let mut secret_generator = SecretRandomGenerator::<ActivatedRandomGenerator>::new(seeder.seed());
    let mut encryption_generator = EncryptionRandomGenerator::<ActivatedRandomGenerator>::new(seeder.seed(), seeder);
    let delta = 1u64 << 58;
    if n < n0 as u64  || m < m0 || m0 > n0 {
        panic!("❗️n0 <= n, m0 <= m, and m0 <= n0 are required");
    }
    if (n0 as u64 + N as u64 / 2) > n {
        panic!("❗️n0 + N/2 <= n is required");
    }
    if let Some(s) = solution {
        if s >= 0 && (s as u64 + m0 as u64 > n0 as u64) {
            panic!("❗️solution + m0 <= n0 is required");
        }
    }
    println!("log_n = {}", log_n);
    println!("n     = {}", n);
    println!("n0    = {}", n0);
    println!("m     = {}", m);
    println!("m0    = {}", m0);
    println!("N     = {}", N);
    println!("sol   = {}", match solution {
        Some(s) => s.to_string(),
        None => solution_val.to_string(),
    });
    println!("Seed  = {}", seed);
    println!("n_nibbles = {}", n_nibbles);
    println!("n_rems    = {}", n_rems);
    let mut T = vec![5u64; n as usize];
    for i in 0..n0 {
        T[i] = rng.random_range(1..=4);
    }
    let mut pattern = vec![54u64; m];
    if solution == Some(-1) {
        for i in 0..m0 {
            pattern[i] = rng.random_range(1..=4);
        }
    } else {
        for i in 0..m0 {
            pattern[i] = T[solution_val as usize + i];
        }
    }
    let mut T_shifted = vec![5u64; n as usize];
    for i in 0..n as usize {
        if i >= N/2 {
            let idx = i - N/2;
            if idx < n as usize {
                T_shifted[i] = T[idx];
            } else {
                T_shifted[i] = 5;
            }
        }
    }
    let symbol = |v| match v {
        1 => "A",
        2 => "C",
        3 => "G",
        4 => "T",
        5 => "$",
        54 => "*",
        _ => "?",
    };
    println!("T         : {}", (0..30).map(|i| symbol(T[i])).collect::<Vec<_>>().join(""));
    println!("T_shifted : {}", (0..30).map(|i| symbol(T_shifted[i])).collect::<Vec<_>>().join(""));
    println!("P         : {}", (0..m.min(30)).map(|i| symbol(pattern[i])).collect::<Vec<_>>().join(""));
    let sa = {
        let n = T.len();
        let mut sa: Vec<usize> = (0..n).collect();
        let mut rank: Vec<usize> = T.iter().map(|&v| v as usize).collect();
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
    let lut_parts = {
        let n = sa.len();
        let mut parts = vec![vec![0u64; n]; n_nibbles];
        for i in 0..n {
            let v = sa[i] as u64;
            for part in 0..n_nibbles {
                parts[part][i] = ((v >> (part * 4)) & 0xF) as u64;
            }
        }
        parts
    };
    println!("========== Key Generation ==========");
    let start_keygen = std::time::Instant::now();   
    let (
        lwe_secret_key,
        glwe_secret_key,
        lwe_secret_key_after_ks,
        bsk,
        ksk,
    ) = keygen_pbs(
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
        auto_fft_type,
        &glwe_secret_key,
        glwe_modular_std_dev,
        &mut encryption_generator,
    );
    let ss_key = generate_scheme_switching_key(
        &glwe_secret_key,
        ss_base_log,
        ss_level,
        glwe_modular_std_dev,
        ciphertext_modulus,
        &mut encryption_generator,
    );
    let mut ct_text = GlweCiphertextList::new(
        0u64,
        glwe_size,
        polynomial_size,
        GlweCiphertextCount((n / N as u64) as usize), 
        ciphertext_modulus,
    );
    let mut ct_text_shifted = GlweCiphertextList::new(
        0u64,
        glwe_size,
        polynomial_size,
        GlweCiphertextCount((n / N as u64) as usize), 
        ciphertext_modulus,
    );
    let mut ct_sa: Vec<GlweCiphertextList<Vec<u64>>> = (0..n_nibbles)
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
    for (glwe_idx, mut glwe) in ct_text.iter_mut().enumerate() {
        let mut plaintext_poly = vec![0u64; polynomial_size.0];
        for (i, val) in plaintext_poly.iter_mut().enumerate() {
            *val = (T[glwe_idx * N + i] as u64) * delta;
        }
        encrypt_glwe_ciphertext(
            &glwe_secret_key,
            &mut glwe,
            &PlaintextList::from_container(plaintext_poly),
            glwe_modular_std_dev,
            &mut encryption_generator,
        );
    }
    for (glwe_idx, mut glwe) in ct_text_shifted.iter_mut().enumerate() {
        let mut plaintext_poly = vec![0u64; polynomial_size.0];
        for (i, val) in plaintext_poly.iter_mut().enumerate() {
            *val = (T_shifted[glwe_idx * N + i] as u64) * delta;
        }
        encrypt_glwe_ciphertext(
            &glwe_secret_key,
            &mut glwe,
            &PlaintextList::from_container(plaintext_poly),
            glwe_modular_std_dev,
            &mut encryption_generator,
        );
    }
    for part in 0..n_nibbles {
        for (glwe_idx, mut cur_glwe) in ct_sa[part].iter_mut().enumerate() {
            let mut plaintext_poly = vec![0u64; polynomial_size.0];
            for (i, val) in plaintext_poly.iter_mut().enumerate() {
                *val = (lut_parts[part][glwe_idx * N + i] as u64) * delta;
            }
            encrypt_glwe_ciphertext(
                &glwe_secret_key,
                &mut cur_glwe,
                &PlaintextList::from_container(plaintext_poly),
                glwe_modular_std_dev,
                &mut encryption_generator,
            );
        }
    }
    let ct_pattern: Vec<LweCiphertext<Vec<u64>>> = pattern.iter().map(|&v| {
        let mut ct = LweCiphertext::new(
            0u64,
            lwe_secret_key.lwe_dimension().to_lwe_size(),
            ciphertext_modulus,
        );
        encrypt_lwe_ciphertext(
            &lwe_secret_key,
            &mut ct,
            Plaintext(v * delta),
            glwe_modular_std_dev,
            &mut encryption_generator,
        );
        ct
    }).collect();
    let keys = FheKeys {
        lwe_secret_key: lwe_secret_key.clone(),
        glwe_secret_key,
        lwe_secret_key_after_ks,
        bsk,
        ksk,
        auto_keys,
        ss_key,
    };
    let fhe_params = FheParams {
        log_n,
        n,
        m,
        n_nibbles,
        n_rems,
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
    println!("🔑 {:6} [ms] in key generation", elapsed_keygen.as_millis());

    /* Server Side */
    let start_find = std::time::Instant::now();   
    println!("========== FHE find ==========");
    let sol = fhe_find(ct_text, ct_text_shifted, ct_sa, ct_pattern, keys, fhe_params);
    let elapsed_find = start_find.elapsed();
    println!("🔍 {:6} [ms] in FHE find", elapsed_find.as_millis());

    /* Client Side (we omit the part of decryption request of sol + r to the Owner, where r is a random value)  */
    let mut sol_recon: u16 = 0;
    for (i, ct) in sol.iter().enumerate() {
        let decrypted = decrypt_lwe_ciphertext(&lwe_secret_key, ct).0;
        let part = ((decrypted + (delta / 2)) / delta) as u16 & 0xF;
        sol_recon |= part << (4 * i);
    }
    let sol_pos = sol_recon as usize;
    if sol_pos + m0 <= T.len() {
        let t_substring: Vec<u64> = T[sol_pos..sol_pos + m0].to_vec();
        let match_found = t_substring == pattern[..m0].to_vec();
        if match_found {
            println!("✅ Match: T[{}:{}] == P", sol_pos, sol_pos + m0);
        } else {
            println!("❌ No Match (error value = {})", sol_recon);
        }
    } else {
        println!("❌ No Match (error value = {})", sol_recon);
    }
}