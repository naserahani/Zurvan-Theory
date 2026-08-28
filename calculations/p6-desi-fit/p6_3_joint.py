#!/usr/bin/env python3
"""P6.3 part 2: calibration, validation, and the joint fit."""

from pathlib import Path

import numpy as np
from scipy.optimize import minimize, brentq
from scipy.interpolate import interp1d

from p6_3_cmb import (E_lcdm_factory, E_hde_factory, chi2_cmb, chi2_bao,
                      r_s, z_star, z_drag, CMB_D, CMB_ICOV, OM_R,
                      chi_comoving, C_KMS, C_OVER_H100)
from p6_2_pantheon import chi2_sn

# ---- calibration: the z_* and z_drag fitting formulas carry known offsets.
# We remove them by one multiplicative factor per observable, fixed at the
# Planck LCDM best fit.  The factors are properties of pre-recombination
# physics, which this model does not modify, so they carry over unchanged.
OM_P, H_P, OB_P = 0.3153, 0.6736, 0.02237
_E = E_lcdm_factory(OM_P, H_P)
_c2, _R, _lA = chi2_cmb(_E, OM_P, H_P, OB_P)
F_R = CMB_D[0] / _R
F_LA = CMB_D[1] / _lA
F_RD = 147.09 / r_s(z_drag(OM_P * H_P**2, OB_P), _E, H_P, OB_P)

print("calibration factors (fixed at Planck LCDM best fit):")
print(f"   f_R  = {F_R:.6f}   f_lA = {F_LA:.6f}   f_rd = {F_RD:.6f}")


def observables(Efun, Om0, h, ob):
    om = Om0 * h**2
    zst = z_star(om, ob)
    DM = chi_comoving(zst, Efun, h)
    R = F_R * np.sqrt(Om0) * h * 100.0 * DM / C_KMS
    lA = F_LA * np.pi * DM / r_s(zst, Efun, h, ob)
    rd = F_RD * r_s(z_drag(om, ob), Efun, h, ob)
    return R, lA, rd


def chi2_all(Om0, h, ob, model, want=("cmb", "bao", "sn")):
    if not (0.05 < Om0 < 0.8 and 0.40 < h < 1.0 and 0.015 < ob < 0.030):
        return 1e8
    try:
        E = E_lcdm_factory(Om0, h) if model == "lcdm" else E_hde_factory(Om0, h)
        R, lA, rd = observables(E, Om0, h, ob)
    except Exception:
        return 1e8
    tot = 0.0
    if "cmb" in want:
        d = np.array([R, lA, ob]) - CMB_D
        tot += float(d @ CMB_ICOV @ d)
    if "bao" in want:
        tot += chi2_bao(E, h, rd)
    if "sn" in want:
        tot += chi2_sn(E)
    return tot


def fit(model, want, x0=(0.31, 0.67, 0.0224)):
    r = minimize(lambda p: chi2_all(p[0], p[1], p[2], model, want), x0,
                 method="Nelder-Mead",
                 options=dict(xatol=1e-6, fatol=1e-7, maxiter=4000))
    return r


if __name__ == "__main__":
    print("\n" + "=" * 72)
    print("  VALIDATION: CMB distance priors alone, flat LCDM")
    print("=" * 72)
    r = fit("lcdm", ("cmb",))
    Om, h, ob = r.x
    print(f"   Omega_m0 = {Om:.4f}   H0 = {100*h:.2f}   omega_b = {ob:.5f}")
    print(f"   omega_m  = {Om*h**2:.5f}      chi2 = {r.fun:.3f} / 3")
    print(f"   [Planck 2018 LCDM: Omega_m = 0.3153, H0 = 67.36, "
          f"omega_m = 0.1430]")

    print("\n" + "=" * 72)
    print("  JOINT FIT: CMB + DESI DR2 BAO + Pantheon+")
    print("=" * 72)
    out = {}
    for model, lab in (("lcdm", "flat LCDM"), ("hde", "HDE c_H=1")):
        r = fit(model, ("cmb", "bao", "sn"))
        Om, h, ob = r.x
        out[lab] = (Om, h, ob, r.fun)
        cmb_only = chi2_all(Om, h, ob, model, ("cmb",))
        bao_only = chi2_all(Om, h, ob, model, ("bao",))
        sn_only = chi2_all(Om, h, ob, model, ("sn",))
        print(f"\n{lab}")
        print(f"   Omega_m0 = {Om:.4f}   H0 = {100*h:.2f}   "
              f"omega_b = {ob:.5f}   omega_m = {Om*h**2:.5f}")
        print(f"   chi2 total = {r.fun:.2f}   "
              f"[CMB {cmb_only:.2f} | BAO {bao_only:.2f} | SN {sn_only:.2f}]")

    dl = out["HDE c_H=1"][3] - out["flat LCDM"][3]
    print(f"\n   Delta chi2 (HDE - LCDM) = {dl:+.2f}   "
          f"[same parameter count: 3]")
    np.save(Path(__file__).with_name("p6_3_out.npy"), out, allow_pickle=True)
