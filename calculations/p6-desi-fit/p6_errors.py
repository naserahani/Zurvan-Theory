#!/usr/bin/env python3
"""P6 part 2: uncertainties, c_H constraint, residuals, figure."""

import numpy as np
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import minimize, brentq

from p6_fit import (hde_E, lcdm_E, chi2, observables, BGS, ANISO,
                    C_OVER_H100, Z_MAX)

N_DATA, DOF = 13, 11


def chi2_hde(Om0, hrd, cH=1.0):
    try:
        E, _ = hde_E(Om0, cH)
    except Exception:
        return 1e6
    return chi2(E, hrd)


def chi2_lcdm(Om0, hrd):
    return chi2(lcdm_E(Om0), hrd)


# ------------------------------------------------- best fits
r_l = minimize(lambda p: chi2_lcdm(*p), [0.30, 101.0], method="Nelder-Mead",
               options=dict(xatol=1e-7, fatol=1e-9, maxiter=5000))
r_h = minimize(lambda p: chi2_hde(p[0], p[1], 1.0), [0.28, 100.0],
               method="Nelder-Mead",
               options=dict(xatol=1e-7, fatol=1e-9, maxiter=5000))

Om_l, hrd_l, c2_l = r_l.x[0], r_l.x[1], r_l.fun
Om_h, hrd_h, c2_h = r_h.x[0], r_h.x[1], r_h.fun


def profile_sigma(f, best, lo, hi):
    """1-sigma from Delta chi2 = 1 on the profile likelihood."""
    g = lambda v: f(v) - f(best) - 1.0
    try:
        left = brentq(g, lo, best, xtol=1e-6)
    except ValueError:
        left = np.nan
    try:
        right = brentq(g, best, hi, xtol=1e-6)
    except ValueError:
        right = np.nan
    return best - left, right - best


# profile over hrd for Om0
def prof_Om_h(Om0):
    r = minimize(lambda h: chi2_hde(Om0, h[0], 1.0), [hrd_h],
                 method="Nelder-Mead", options=dict(xatol=1e-6, fatol=1e-9))
    return r.fun

def prof_Om_l(Om0):
    r = minimize(lambda h: chi2_lcdm(Om0, h[0]), [hrd_l],
                 method="Nelder-Mead", options=dict(xatol=1e-6, fatol=1e-9))
    return r.fun

sm_l = profile_sigma(prof_Om_l, Om_l, 0.20, 0.42)
sm_h = profile_sigma(prof_Om_h, Om_h, 0.18, 0.40)


def prof_hrd_h(hrd):
    r = minimize(lambda o: chi2_hde(o[0], hrd, 1.0), [Om_h],
                 method="Nelder-Mead", options=dict(xatol=1e-6, fatol=1e-9))
    return r.fun

def prof_hrd_l(hrd):
    r = minimize(lambda o: chi2_lcdm(o[0], hrd), [Om_l],
                 method="Nelder-Mead", options=dict(xatol=1e-6, fatol=1e-9))
    return r.fun

sh_l = profile_sigma(prof_hrd_l, hrd_l, 95.0, 108.0)
sh_h = profile_sigma(prof_hrd_h, hrd_h, 93.0, 107.0)

# ------------------------------------------------- c_H constraint
def prof_cH(cH):
    r = minimize(lambda p: chi2_hde(p[0], p[1], cH), [Om_h, hrd_h],
                 method="Nelder-Mead", options=dict(xatol=1e-6, fatol=1e-9,
                                                    maxiter=3000))
    return r.fun

r_c = minimize(lambda p: chi2_hde(p[0], p[1], p[2]), [Om_h, hrd_h, 1.0],
               method="Nelder-Mead",
               options=dict(xatol=1e-7, fatol=1e-9, maxiter=8000))
cH_best, c2_c = r_c.x[2], r_c.fun
s_c = profile_sigma(prof_cH, cH_best, 0.55, 1.8)

dev_cH = abs(1.0 - cH_best) / (s_c[1] if cH_best < 1 else s_c[0])

# ------------------------------------------------- w(z)
_, OLf = hde_E(Om_h, 1.0)
w = lambda z: -1.0/3.0 - (2.0/3.0)*np.sqrt(float(OLf(z)))
w0 = w(0.0)
zp = 0.5
wa_eff = (w(zp) - w0) / (zp/(1+zp))

print("=" * 70)
print(" P6 RESULTS")
print("=" * 70)
print(f"\nflat LCDM      Om0  = {Om_l:.4f} (-{sm_l[0]:.4f} +{sm_l[1]:.4f})")
print(f"               h*rd = {hrd_l:.2f} (-{sh_l[0]:.2f} +{sh_l[1]:.2f}) Mpc")
print(f"               chi2 = {c2_l:.3f} / {DOF} = {c2_l/DOF:.3f}")
print(f"   [DESI DR2 official: Om0 = 0.2975 +/- 0.0086,  h*rd = 101.54 +/- 0.73]")

print(f"\nHDE c_H = 1    Om0  = {Om_h:.4f} (-{sm_h[0]:.4f} +{sm_h[1]:.4f})")
print(f"               h*rd = {hrd_h:.2f} (-{sh_h[0]:.2f} +{sh_h[1]:.2f}) Mpc")
print(f"               chi2 = {c2_h:.3f} / {DOF} = {c2_h/DOF:.3f}")
print(f"               w0   = {w0:.4f}   w_a(eff) = {wa_eff:+.4f}")

print(f"\nDelta chi2 (HDE_c1 - LCDM) = {c2_h - c2_l:+.3f}  [same param count]")

print(f"\nc_H free:      c_H  = {cH_best:.4f} (-{s_c[0]:.4f} +{s_c[1]:.4f})")
print(f"               chi2 = {c2_c:.3f}")
print(f"               c_H = 1 lies {dev_cH:.2f} sigma from best fit")
print(f"               Delta chi2 for fixing c_H = 1: {c2_h - c2_c:.3f}")

# ------------------------------------------------- residuals
print("\nResiduals (obs - model)/sigma, c_H = 1 fit:")
zl = [BGS["z"]] + [t["z"] for t in ANISO]
Eh, _ = hde_E(Om_h, 1.0)
oh = observables(Eh, hrd_h, zl)
ol = observables(lcdm_E(Om_l), hrd_l, zl)
print(f"  {'tracer':<11}{'z':>6}{'obs':>10}{'HDE':>9}{'pull':>7}{'LCDM':>9}{'pull':>7}")
d = oh[BGS['z']][2]; e = ol[BGS['z']][2]
print(f"  {'BGS  DV':<11}{BGS['z']:>6.3f}{BGS['DV']:>10.3f}{d:>9.3f}"
      f"{(BGS['DV']-d)/BGS['sV']:>7.2f}{e:>9.3f}{(BGS['DV']-e)/BGS['sV']:>7.2f}")
for t in ANISO:
    for k, lbl in ((0, "DM"), (1, "DH")):
        mh, ml = oh[t["z"]][k], ol[t["z"]][k]
        obs_, s = (t["DM"], t["sM"]) if k == 0 else (t["DH"], t["sH"])
        print(f"  {t['name']+'  '+lbl:<11}{t['z']:>6.3f}{obs_:>10.3f}{mh:>9.3f}"
              f"{(obs_-mh)/s:>7.2f}{ml:>9.3f}{(obs_-ml)/s:>7.2f}")

# ------------------------------------------------- figure
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
zz = np.linspace(0.01, Z_MAX, 300)

Kh = C_OVER_H100 / hrd_h
Kl = C_OVER_H100 / hrd_l
ohz = observables(Eh, hrd_h, zz)
olz = observables(lcdm_E(Om_l), hrd_l, zz)
dm_h = np.array([ohz[z][0] for z in zz]); dm_l = np.array([olz[z][0] for z in zz])
dh_h = np.array([ohz[z][1] for z in zz]); dh_l = np.array([olz[z][1] for z in zz])

ax[0].plot(zz, dm_h/dm_l - 1, lw=2, color="#b5533c", label=r"HDE $c_H=1$")
ax[0].axhline(0, color="0.35", lw=1.5, ls="--", label=r"flat $\Lambda$CDM")
for t in ANISO:
    m = np.interp(t["z"], zz, dm_l)
    ax[0].errorbar(t["z"], t["DM"]/m - 1, yerr=t["sM"]/m, fmt="o",
                   ms=4.5, color="#2a3f5f", capsize=2.5, lw=1.2)
ax[0].set_xlabel("z"); ax[0].set_ylabel(r"$(D_M/r_d)\,/\,\Lambda$CDM $-\,1$")
ax[0].set_title("Transverse distance", fontsize=10)
ax[0].legend(frameon=False, fontsize=9)

ax[1].plot(zz, [w(z) for z in zz], lw=2, color="#b5533c",
           label=r"HDE $c_H=1$")
ax[1].axhline(-1, color="0.35", lw=1.5, ls="--", label=r"$\Lambda$CDM")
ax[1].plot(zz, -0.42 - 1.75*zz/(1+zz), lw=1.6, color="#4a7c59", ls=":",
           label=r"DESI CPL ($w_0{=}-0.42$, $w_a{=}-1.75$)")
ax[1].set_xlabel("z"); ax[1].set_ylabel(r"$w(z)$")
ax[1].set_title("Equation of state — note opposite sign of $w_a$", fontsize=10)
ax[1].legend(frameon=False, fontsize=8, loc="lower left")
ax[1].set_ylim(-2.1, -0.35)

for a in ax:
    a.grid(alpha=0.18, lw=0.6)
    for s in ("top", "right"):
        a.spines[s].set_visible(False)
plt.tight_layout()
plt.savefig(Path(__file__).with_name("p6_hde_vs_desi.png"), dpi=170)
print("\nfigure saved")
