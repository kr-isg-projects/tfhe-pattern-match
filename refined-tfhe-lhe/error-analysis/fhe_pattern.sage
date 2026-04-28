# This file (var.sage) is a modified version of:
# https://github.com/KAIST-CryptLab/refined-tfhe-lhe/blob/main/error_analysis/integer_input_lhe.sage
# Original file source: KAIST-CryptLab/refined-tfhe-lhe
# For review or reference purposes only

load("var.sage")

stddev_653 = 6.772274609789095e-05 # 130.6 bit security
stddev_636 = 9.25119974676756e-5 # 130.7 bit security
stddev_768 = 8.763872947670246e-06 # 130.1 bit security
stddev_769 = 8.763872947670246e-06 # 130.1 bit security
stddev_808 = 2.1124945159091033e-06 # 128 bit security
stddev_873 = 1.3962609252411138e-06 # 130.1 bit security
stddev_953 = 3.4906523131027844e-07 # 130.2 bit security
stddev_2048 = 9.25119974676756e-16 # 130.7 bit security
stddev_4096 = 2.168404344971009e-19 # 218.5 bit security

FHE_KS_PBS = (
    "FHE_KS_PBS",
    (808, stddev_808^2), # (LWE dim, LWE var),
    (2^11, 3), # (PBS base, PBS level) (2^15, 2) 
    (2^12, 4, 2^40), # (Tr base, Tr level, split fft)
    (2^10, 4), # (SS base, SS level),
    (2^5, 4), # (CBS base, CBS level)
    (2^2, 8), # (KS base, KS level)
    6, # log_modulus
    (4096, stddev_4096^2), # (N, Var_GLWE)
)
q = 2^64
k = 1
param = FHE_KS_PBS
name = param[0]
(n, Var_LWE) = param[1]
(B_pbs, l_pbs) = param[2]
(B_tr, l_tr, b_tr) = param[3]
(B_ss, l_ss) = param[4]
(B_cbs, l_cbs) = param[5]
(B_ksk, l_ksk) = param[6]
log_modulus = param[7]
(N, Var_GLWE) = param[8]
print(f"========================= {name} =========================")
print(f"n: {n}, N: {N}, k: {k}")
print(f"B_ksk: 2^{log(B_ksk, 2)}, l_ksk: {l_ksk}")
print(f"B_pbs: 2^{log(B_pbs, 2)}, l_pbs: {l_pbs}")
print(f"B_tr: 2^{log(B_tr, 2)}, l_tr: {l_tr}, b_tr: 2^{log(b_tr, 2)}")
print(f"B_ss: 2^{log(B_ss, 2)}, l_ss: {l_ss}")
print(f"B_cbs: 2^{log(B_cbs, 2)}, l_cbs: {l_cbs}")
print(f"Var_LWE: 2^{log(q^2 * Var_LWE, 2).n():.4f}") # used only after KS and before PBS
print(f"Var_GLWE: 2^{log(q^2 * Var_GLWE, 2).n():.4f}") # = Var_LWE except after KS

print("-------------------------------------------------------------------------")
print("Noise Variance of CBS & XP in [WHS25]")
Var_pbs = get_var_pbs(N, k, n, q, Var_GLWE, B_pbs, l_pbs)
Var_fft_pbs = get_var_fft_pbs(N, k, n, B_pbs, l_pbs)
Var_pbs_tot = Var_pbs + Var_fft_pbs
Var_ks = get_var_lwe_ks(N, k, q, Var_LWE, B_ksk, l_ksk) # N0 to N 

Var_tr = get_var_tr(N, k, q, Var_GLWE, B_tr, l_tr)
Var_auto_gadget = get_var_glwe_ks_gadget(N, k, q, B_tr, l_tr)
Var_auto_key = get_var_glwe_ks_key(N, k, q, Var_GLWE, B_tr, l_tr)
Var_fft_tr = get_var_fft_tr(N, k, B_tr, l_tr, b_tr)
Var_tr_tot = Var_tr + Var_fft_tr

Var_split_fft_tr_upper = get_var_fft_glwe_ks(N, k, B_tr, l_tr, q / b_tr)
_, fp_split_fft = get_fp_split_fft_glwe_ks(N, k, q, B_tr, l_tr, b_tr)
log_fp_split_fft_tr = log(fp_split_fft, 2).n(10000)

Var_ss = get_var_ss(N, k, q, q^2 * Var_GLWE, B_ss, l_ss)
Var_fft_ss = get_var_fft_ext_prod(N, k, q, B_ss, l_ss)
Var_ss_tot = Var_ss + Var_fft_ss

Var_cbs = Var_pbs_tot + Var_ss_tot + (N/2) * (Var_tr_tot )
Var_cbs_additive = Var_pbs_tot + Var_ss_tot
Var_cbs_amp = (N/2) * (Var_tr_tot )

Var_ext = get_var_ext_prod(N, k, q, Var_cbs, B_cbs, l_cbs) + get_var_fft_ext_prod(N, k, q, B_cbs, l_cbs)
Var_add_gadget = get_var_ext_prod_gadget(N, k, q, B_cbs, l_cbs)
Var_add_inc = get_var_ext_prod_inc(N, k, Var_cbs, B_cbs, l_cbs)
Var_fft_ext = get_var_fft_ext_prod(N, k, q, B_cbs, l_cbs)

print(f"Var_ks: 2^{log(Var_ks, 2).n():.4f}")
print(f"Var_pbs_tot: 2^{log(Var_pbs_tot, 2).n():.4f}")
print(f"Var_tr_tot : 2^{log(Var_tr_tot, 2).n():.4f}")
print(f"  - F.P. of split fft: 2^{log_fp_split_fft_tr:.4f} (stddev_upper: 2^{log(Var_split_fft_tr_upper, 2).n() / 2:.4f})")
print(f"Var_ss_tot: 2^{log(Var_ss_tot, 2).n():.4f}")
print(f"Var_cbs: 2^{log(Var_cbs, 2).n():.4f}")
print(f"Var_ext: 2^{log(Var_ext, 2).n():.4f}")
Var_LWE = q^2 * Var_LWE
log_fp_thrs = -40

print("-------------------------------------------------------------------------")
print("Noise Analysis in our circuit (LHE mode)")
theta = log(l_cbs, 2) # PBSmanyLUT with 2^theta PBS
delta_in = 2^(64 - 1) # delta in for PBSmanyLUT
print("theta (we have 2^theta) for PBSmanyLUT:", theta) # LUTmanyPBS with 2^theta PBS
print(f"Delta_in for PBSmanyLUT: 2^{log(delta_in, 2)}")
_, min_fp = get_min_fp_pbs(n, q, N, theta, delta_in)
log_min_fp = log(min_fp, 2).n(1000)
print(f"Min F.P.: 2^{log_min_fp:.4f}")
if log_min_fp > log_fp_thrs:
    print(f"  - Var_thrs_{-log_fp_thrs:.0f}: impossible")
else:
    log_Var_thrs = find_var_thrs(n, q, N, theta, delta_in, log_fp_thrs)
    Var_thrs = 2^log_Var_thrs
    Var_in = Var_pbs_tot
    Var_circ = 2^(2*(log_modulus - 1)) * 2 * Var_in + Var_ks
    _,fp = get_fp_pbs(n, q, N, theta, delta_in, Var_in=Var_circ)
    log_fp = log(fp, 2).n(1000)
    print(f"Var_thrs: 2^{log(Var_thrs.n(), 2)}")
    print(f"Var_circ: 2^{log(Var_circ.n(), 2)}")
    print(f"F.P. for PBSmanyLUT: 2^{log_fp:.4f}")

print("-------------------------------------------------------------------------")
print("Noise Analysis in our circuit (FHE mode)")
theta = 0 # normal PBS
delta_in = 2^(64 - 6) # delta in for PBSmanyLUT
print("theta (we have 2^theta) for PBSmanyLUT:", theta) # LUTmanyPBS with 2^theta PBS
print(f"Delta_in for PBS: 2^{log(delta_in, 2)}")
_, min_fp = get_min_fp_pbs(n, q, N, theta=theta, delta_in=delta_in)
log_min_fp = log(min_fp, 2).n(1000)
print(f"Min F.P.: 2^{log_min_fp:.4f}")
if log_min_fp > log_fp_thrs:
    print(f"  - Var_thrs_{-log_fp_thrs:.0f}: impossible")
else:
    log_Var_thrs = find_var_thrs(n, q, N, theta, delta_in, log_fp_thrs)
    Var_thrs = 2^log_Var_thrs
    Var_in = Var_pbs_tot
    Var_circ = 2*Var_in + Var_tr_tot + Var_ext + Var_ks
    _,fp = get_fp_pbs(n, q, N, theta, delta_in, Var_in=Var_circ)
    log_fp = log(fp, 2).n(1000)
    print(f"Var_thrs:2^{log(Var_thrs.n(), 2)}")
    print(f"Var_circ:2^{log(Var_circ.n(), 2)}")
    print(f"F.P. for PBS: 2^{log_fp:.4f}")