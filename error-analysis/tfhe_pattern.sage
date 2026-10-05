# This file is a modified version of:
# https://github.com/KAIST-CryptLab/refined-tfhe-lhe/blob/main/error_analysis/integer_input_lhe.sage
# Original file source: KAIST-CryptLab/refined-tfhe-lhe

load("var.sage")

stddev_808  = 3.588664896844116206e-05  # 128 bit security
stddev_920  = 6.843470867416183799e-06  # 128 bit security
stddev_2048 = 3.54e-13 # 128 bit security
stddev_4096 = 3.014727391206267763e-23 # 128 bit security

FHE_KS_PBS = (
    "TFHE_RS_2_2_WOPBS",
    (920, stddev_920^2),
    (2^11, 3),         # PBS
    (2^12, 2, 2^40),   # HomTrace / GLWE KS (our circuit)
    (2^2, 7),          # CBS / External Product
    (2^2, 7),          # LWE KS
    (2^12, 2),         # PFPKS
    4,                 # log2 plaintext modulus; p = 16
    (2048, stddev_2048^2),
)

q = 2^64
k = 1
param = FHE_KS_PBS
log_fp_thrs = -128

name = param[0]
(n, Var_LWE) = param[1]
(B_pbs, l_pbs) = param[2]
(B_tr, l_tr, b_tr) = param[3]
(B_cbs, l_cbs) = param[4]
(B_ksk, l_ksk) = param[5]
(B_pfpks, l_pfpks) = param[6]
log_modulus = param[7]
(N, Var_GLWE) = param[8]

print(f"========================= {name} =========================")
print(f"n: {n}, N: {N}, k: {k}")
print(f"B_ksk: 2^{log(B_ksk, 2)}, l_ksk: {l_ksk}")
print(f"B_pbs: 2^{log(B_pbs, 2)}, l_pbs: {l_pbs}")
print(f"B_tr: 2^{log(B_tr, 2)}, l_tr: {l_tr}, b_tr: {b_tr}")
print(f"B_cbs: 2^{log(B_cbs, 2)}, l_cbs: {l_cbs}")
print(f"B_pfpks: 2^{log(B_pfpks, 2)}, l_pfpks: {l_pfpks}")
print(f"Var_LWE: 2^{log(q^2 * Var_LWE, 2).n():.4f}")
print(f"Var_GLWE: 2^{log(q^2 * Var_GLWE, 2).n():.4f}")

print("-------------------------------------------------------------------------")
print("Noise Variance of native CBS & External Product")

Var_pbs = get_var_pbs(N, k, n, q, Var_GLWE, B_pbs, l_pbs)
Var_fft_pbs = get_var_fft_pbs(N, k, n, B_pbs, l_pbs)
Var_pbs_tot = Var_pbs + Var_fft_pbs
Var_ks = get_var_lwe_ks(N, k, q, Var_LWE, B_ksk, l_ksk)
Var_cbs = get_var_cbs(
    N, k, n, q, Var_GLWE,
    B_pbs, l_pbs,
    B_pfpks, l_pfpks
)
Var_ext = get_var_ext_prod(
    N, k, q, Var_cbs, B_cbs, l_cbs
) + get_var_fft_ext_prod(
    N, k, q, B_cbs, l_cbs
)

print(f"Var_ks: 2^{log(Var_ks, 2).n():.4f}")
print(f"Var_pbs_tot: 2^{log(Var_pbs_tot, 2).n():.4f}")
print(f"Var_cbs: 2^{log(Var_cbs, 2).n():.4f}")
print(f"Var_ext: 2^{log(Var_ext, 2).n():.4f}")

# -------- Noise variance of the longest path before PBS -------- #
# LWE noise after PBS
Var_lwe = Var_pbs_tot

# HomTrace / RLWE KS
Var_tr = get_var_tr(N, k, q, Var_GLWE, B_tr, l_tr)
Var_fft_tr = get_var_fft_tr(N, k, B_tr, l_tr, b_tr)
Var_tr_tot = Var_tr + Var_fft_tr

# CMux / External Product
Var_xp = get_var_ext_prod(
    N, k, q, Var_cbs, B_cbs, l_cbs
) + get_var_fft_ext_prod(
    N, k, q, B_cbs, l_cbs
)

# LWE KS before PBS
Var_ks = get_var_lwe_ks(
    N, k, q, Var_LWE, B_ksk, l_ksk
)

# 2*V_lwe + V_tr + V_xp + V_ks
Var_path_pbs = 2*Var_lwe + Var_tr_tot + Var_xp + Var_ks

print(f"Var_lwe:     2^{log(Var_lwe, 2).n():.4f}")
print(f"Var_tr_tot:  2^{log(Var_tr_tot, 2).n():.4f}")
print(f"Var_xp:      2^{log(Var_xp, 2).n():.4f}")

Var_ks_gadget = get_var_lwe_ks_gadget(N, k, q, B_ksk, l_ksk)
Var_ks_key = get_var_lwe_ks_key(N, k, q, Var_LWE, B_ksk, l_ksk)
print(f"Var_ks_gadget: 2^{log(Var_ks_gadget, 2).n():.4f}")
print(f"Var_ks_key:    2^{log(Var_ks_key, 2).n():.4f}")
print(f"Var_ks:        2^{log(Var_ks, 2).n():.4f}")
print(f"Var_path_pbs: 2^{log(Var_path_pbs, 2).n():.4f}")


# -------- Longest noise paths -------- #
p = 2^log_modulus

# Path 1: l' + r' -> HomTrace -> CMux -> SE -> LWE KS -> PBS
Var_path_pbs = 2*Var_lwe + Var_tr_tot + Var_xp + Var_ks

# Path 2: two LWE additions -> shift by at most p/2 -> LWE KS -> native CBS
Var_path_cbs = (p^2/2)*Var_lwe + Var_ks

# Path 3: CMux tree -> BR
# Set this to the maximum supported padded text capacity
T_size = 2^20
cmux_depth = max(0, ceil(log(T_size/N, 2)))
Var_path_cmuxtree = cmux_depth * Var_xp

Var_circ = max(
    Var_path_pbs,
    Var_path_cbs,
    Var_path_cmuxtree
)

print("-------------------------------------------------------------------------")
print("Longest Noise Paths")
print(f"CMux tree depth: {cmux_depth}")
print(f"Var_path_pbs:      2^{log(Var_path_pbs, 2).n():.4f}")
print(f"Var_path_cbs:      2^{log(Var_path_cbs, 2).n():.4f}")

if cmux_depth > 0:
    print(f"Var_path_cmuxtree: 2^{log(Var_path_cmuxtree, 2).n():.4f}")
else:
    print("Var_path_cmuxtree: 0")
print(f"Var_circ:          2^{log(Var_circ, 2).n():.4f}")


print("-------------------------------------------------------------------------")
print("Failure Probability Analysis")

theta = 0 # normal PBS
delta_in = 2^(64 - log_modulus)

print(f"Delta_in for PBS: 2^{log(delta_in, 2)}")

_, min_fp = get_min_fp_pbs(n, q, N, theta, delta_in)
log_min_fp = log(min_fp, 2).n(1000)

print(f"Min F.P.: 2^{log_min_fp:.4f}")

if log_min_fp > log_fp_thrs:
    print(f"Var_thrs_{-log_fp_thrs:.0f}: impossible")
else:
    log_Var_thrs = find_var_thrs(
        n, q, N, theta, delta_in, log_fp_thrs
    )
    Var_thrs = 2^log_Var_thrs

    _, fp = get_fp_pbs(
        n, q, N, theta, delta_in, Var_in=Var_circ
    )
    log_fp = log(fp, 2).n(1000)

    print(f"Var_circ:     2^{log(Var_circ, 2).n():.4f}")
    print(f"F.P.:         2^{log_fp:.4f}")

    if Var_circ < Var_thrs:
        print("Overall circuit: OK")
    else:
        print("Overall circuit: NG")
