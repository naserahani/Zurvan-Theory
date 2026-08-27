#!/usr/bin/env python3
"""
P6.2 — Does the c_H = 1 model's preference for a low Omega_m survive the
Pantheon+ supernova sample?

Why supernovae.  DESI BAO constrains only (Omega_m0, h*r_d); a prior on r_d is
orthogonal to that and cannot move Omega_m0 (see p6_1_fast.py).  Supernovae
constrain the *shape* of the distance-redshift relation with no sound horizon
and no LambdaCDM-derived input, so they give an independent handle on Omega_m0.

Method.  Standard SN-only cosmology: fit m_b_corr with the absolute magnitude
(and hence H0) absorbed into a single nuisance offset script_M, analytically
marginalised with a flat prior:

    chi2 = A - B^2/C,   A = r^T C^-1 r,  B = r^T C^-1 1,  C = 1^T C^-1 1

with r = m_obs - 5 log10[(1+z_hel) * Integral_0^zHD dz/E(z)].
Calibrator SNe are excluded and z_HD > 0.01 is imposed.

Data: Pantheon+ (Brout et al. 2022), STAT+SYS covariance, official release.
"""

from pathlib import Path
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import brentq, minimize_scalar

from p6_fit import hde_E, lcdm_E, BGS, ANISO, C_OVER_H100

HERE = Path(__file__).parent / "pantheon"
ZMIN = 0.01

# ------------------------------------------------------------------ load SNe
raw = np.genfromtxt(HERE / "data.txt", names=True, dtype=None, encoding="utf-8")
mask = (raw["zHD"] > ZMIN) & (raw["IS_CALIBRATOR"] == 0)
idx = np.where(mask)[0]

zHD = raw["zHD"][mask]
zHEL = raw["zHEL"][mask]
mb = raw["m_b_corr"][mask]
N = len(zHD)

flat = np.loadtxt(HERE / "cov.txt", skiprows=1)
NFULL = int(open(HERE / "cov.txt").readline())
COV = flat.reshape(NFULL, NFULL)[np.ix_(idx, idx)]
CINV = np.linalg.inv(COV)
ONES = np.ones(N)
C_ONE = float(ONES @ CINV @ ONES)

print(f"Pantheon+: {NFULL} SNe in release, {N} used "
      f"(z > {ZMIN}, calibrators excluded)")
print(f"           z range {zHD.min():.4f} - {zHD.max():.4f}")

# ------------------------------------------------------------------ SN chi2
def chi2_sn(Efun):
    zg = np.linspace(0.0, max(zHD.max(), 2.6) * 1.02, 6000)
    Ev = Efun(zg)
    integ = np.concatenate([[0.0],
                            np.cumsum(np.diff(zg) * 0.5 * (1/Ev[:-1] + 1/Ev[1:]))])
    chi = interp1d(zg, integ, kind="cubic")
    dl_shape = (1.0 + zHEL) * chi(zHD)          # D_L in units of c/H0
    model = 5.0 * np.log10(dl_shape)
    r = mb - model
    A = float(r @ CINV @ r)
    B = float(r @ CINV @ ONES)
    return A - B * B / C_ONE                    # marginalised over script_M


# ------------------------------------------------------------------ BAO chi2
ZS = [BGS["z"]] + [t["z"] for t in ANISO]

def _bao_vectors(Efun):
    zg = np.linspace(0.0, 2.6, 4000)
    Ev = Efun(zg)
    integ = np.concatenate([[0.0],
                            np.cumsum(np.diff(zg) * 0.5 * (1/Ev[:-1] + 1/Ev[1:]))])
    chi = interp1d(zg, integ, kind="cubic")
    d, m, blocks = [], [], []
    z = BGS["z"]; cm = float(chi(z)); E = float(Efun(z))
    d.append(BGS["DV"]); m.append((cm**2 * z / E) ** (1/3))
    blocks.append((1, np.array([[1 / BGS["sV"]**2]])))
    for t in ANISO:
        z = t["z"]; cm = float(chi(z)); E = float(Efun(z))
        d.extend([t["DM"], t["DH"]]); m.extend([cm, 1.0 / E])
        blocks.append((2, t["icov"]))
    return np.array(d), np.array(m), blocks


def chi2_bao(Efun):
    """Minimised analytically over h*r_d (chi2 is quadratic in u = c/(100 h rd))."""
    d, m, blocks = _bao_vectors(Efun)
    num = den = 0.0; i = 0
    for n, ic in blocks:
        dv, mv = d[i:i+n], m[i:i+n]
        num += float(mv @ ic @ dv); den += float(mv @ ic @ mv); i += n
    u = num / den
    r = d - u * m
    out = 0.0; i = 0
    for n, ic in blocks:
        v = r[i:i+n]; out += float(v @ ic @ v); i += n
    return out, C_OVER_H100 / u


# ------------------------------------------------------------------ scan
def scan(model, oms, cH=1.0):
    sn, bao, hrd = [], [], []
    for om in oms:
        E = lcdm_E(om) if model == "lcdm" else hde_E(om, cH)[0]
        sn.append(chi2_sn(E))
        b, x = chi2_bao(E)
        bao.append(b); hrd.append(x)
    return np.array(sn), np.array(bao), np.array(hrd)


def chi2_at(om, model, which, cH=1.0):
    """chi2 for one probe at a single Omega_m0 (no grid)."""
    E = lcdm_E(om) if model == "lcdm" else hde_E(om, cH)[0]
    if which == "sn":
        return chi2_sn(E)
    if which == "bao":
        return chi2_bao(E)[0]
    return chi2_sn(E) + chi2_bao(E)[0]


def minimum(model, which, lo=0.16, hi=0.50, cH=1.0):
    """Continuous minimisation + 1-sigma interval.

    A grid is NOT used here: the true minimum never lands on a grid node, so
    grid output is quantised at the step size (a 0.005 grid shifts the LCDM
    SN/BAO tension from 1.75 to 1.96 sigma).  Bounded Brent finds the minimum
    to 1e-5 in Omega_m0 with far fewer evaluations than a fine grid.
    """
    f = lambda om: chi2_at(om, model, which, cH)
    r = minimize_scalar(f, bounds=(lo, hi), method="bounded",
                        options=dict(xatol=1e-5))
    om0, c0 = r.x, r.fun
    g = lambda om: f(om) - c0 - 1.0
    a = brentq(g, max(lo, om0 - 0.12), om0, xtol=1e-5)
    b = brentq(g, om0, min(hi, om0 + 0.12), xtol=1e-5)
    return om0, om0 - a, b - om0, c0


OMS = np.linspace(0.15, 0.55, 81)      # grid: used ONLY for the figure

print("\n" + "=" * 74)
print("  P6.2 - Pantheon+ supernovae, and combination with DESI DR2 BAO")
print("=" * 74)

store, best = {}, {}
for label, model in (("flat LCDM", "lcdm"), ("HDE c_H=1", "hde")):
    sn, bao, hrd = scan(model, OMS)
    store[label] = dict(sn=sn, bao=bao)
    print(f"\n{label}")
    for which, dof in (("sn", N - 1), ("bao", 11), ("comb", N + 13 - 2)):
        om, ms, ps, c0 = minimum(model, which)
        best[(label, which)] = (om, (ms + ps) / 2, c0)
        print(f"   {which:<9} Omega_m0 = {om:.4f} (-{ms:.4f} +{ps:.4f})"
              f"   chi2 = {c0:.2f} / {dof}")

print("\n" + "-" * 74)
print("Internal tension between SN and BAO, within each model:")
for label in ("flat LCDM", "HDE c_H=1"):
    os_, ss, _ = best[(label, "sn")]
    ob, sb, _ = best[(label, "bao")]
    print(f"   {label:<11} SN {os_:.4f} vs BAO {ob:.4f}"
          f"   -> {abs(os_ - ob) / np.hypot(ss, sb):.2f} sigma")

d = best[("HDE c_H=1", "comb")][2] - best[("flat LCDM", "comb")][2]
print(f"\nCombined Delta chi2 (HDE - LCDM) = {d:+.2f}   [same parameter count]")
# No significance is quoted: the two models are not nested, so Wilks'
# theorem does not apply and Delta chi2 cannot be converted to sigma.
print(f"   (BAO alone: {best[('HDE c_H=1','bao')][2]-best[('flat LCDM','bao')][2]:+.2f};"
      f"  SN alone: {best[('HDE c_H=1','sn')][2]-best[('flat LCDM','sn')][2]:+.2f})")

np.savez(Path(__file__).with_name("p6_2_scan.npz"), oms=OMS,
         lcdm_sn=store["flat LCDM"]["sn"], lcdm_bao=store["flat LCDM"]["bao"],
         hde_sn=store["HDE c_H=1"]["sn"], hde_bao=store["HDE c_H=1"]["bao"])
print("\nscan saved")
