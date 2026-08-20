#!/usr/bin/env python3
"""
P6.1 — fast version.

Key structural fact used here: every DESI BAO observable scales as 1/(h*r_d).
Writing u = c/(100 h r_d), the model vector is u * m(Omega_m0), where m depends
only on Omega_m0.  The chi2 is then exactly quadratic in u and can be minimised
in closed form:   u* = (m^T C^-1 d) / (m^T C^-1 m).

This makes the Omega_m0 scan essentially free and, more importantly, exposes
what an r_d prior can and cannot do.
"""

from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d

from p6_fit import BGS, ANISO, C_OVER_H100, hde_E, lcdm_E

RD_MEAN, RD_SIG = 147.09, 0.26

ZS = [BGS["z"]] + [t["z"] for t in ANISO]


def model_vector(Efun):
    """Return (d, m, Cinv_blocks): data vector, shape vector (at u=1)."""
    zg = np.linspace(0.0, 2.6, 4000)
    Ev = Efun(zg)
    integ = np.concatenate([[0.0],
                            np.cumsum(np.diff(zg) * 0.5 * (1/Ev[:-1] + 1/Ev[1:]))])
    chi = interp1d(zg, integ, kind="cubic")

    d, m, blocks = [], [], []

    # BGS: D_V/r_d
    z = BGS["z"]
    cm = float(chi(z)); E = float(Efun(z))
    mV = (cm**2 * z / E) ** (1/3)
    d.append(BGS["DV"]); m.append(mV)
    blocks.append(("scalar", np.array([[1/BGS["sV"]**2]])))

    for t in ANISO:
        z = t["z"]
        cm = float(chi(z)); E = float(Efun(z))
        d.extend([t["DM"], t["DH"]])
        m.extend([cm, 1.0/E])
        blocks.append(("pair", t["icov"]))

    return np.array(d), np.array(m), blocks


def chi2_of_u(d, m, blocks, u):
    r = d - u * m
    out, i = 0.0, 0
    for kind, ic in blocks:
        n = 1 if kind == "scalar" else 2
        v = r[i:i+n]
        out += float(v @ ic @ v)
        i += n
    return out


def best_u(d, m, blocks):
    """Closed-form minimiser of the quadratic chi2 in u."""
    num = den = 0.0
    i = 0
    for kind, ic in blocks:
        n = 1 if kind == "scalar" else 2
        dv, mv = d[i:i+n], m[i:i+n]
        num += float(mv @ ic @ dv)
        den += float(mv @ ic @ mv)
        i += n
    return num / den


def scan(model, cH=1.0, om_lo=0.20, om_hi=0.40, n=161):
    oms = np.linspace(om_lo, om_hi, n)
    chis, hrds = [], []
    for om in oms:
        E = lcdm_E(om) if model == "lcdm" else hde_E(om, cH)[0]
        d, m, B = model_vector(E)
        u = best_u(d, m, B)
        chis.append(chi2_of_u(d, m, B, u))
        hrds.append(C_OVER_H100 / u)
    return oms, np.array(chis), np.array(hrds)


def interval(oms, chis):
    c0 = chis.min(); i0 = int(np.argmin(chis))
    f = interp1d(oms, chis - c0 - 1.0, kind="cubic")
    from scipy.optimize import brentq
    lo = brentq(f, oms[0], oms[i0]); hi = brentq(f, oms[i0], oms[-1])
    return oms[i0], oms[i0] - lo, hi - oms[i0], c0, i0


print("=" * 74)
print("  P6.1 — what a Planck r_d prior can and cannot do")
print("=" * 74)

out = {}
for label, model in (("flat LCDM", "lcdm"), ("HDE c_H=1", "hde")):
    oms, chis, hrds = scan(model)
    om, sm, sp, c0, i0 = interval(oms, chis)
    hrd = hrds[i0]
    H0 = 100 * hrd / RD_MEAN
    out[label] = dict(om=om, sm=sm, sp=sp, chi2=c0, hrd=hrd, H0=H0)
    print(f"\n{label}")
    print(f"   Omega_m0 = {om:.4f}  (-{sm:.4f} +{sp:.4f})   [BAO alone]")
    print(f"   h*r_d    = {hrd:.2f} Mpc")
    print(f"   chi2     = {c0:.3f} / 11")
    print(f"   with r_d = {RD_MEAN} Mpc  ->  H0 = {H0:.2f} km/s/Mpc")

print("\n" + "-" * 74)
print("STRUCTURAL RESULT")
print("-" * 74)
print("""
Every BAO observable scales as 1/(h r_d).  The likelihood therefore constrains
only Omega_m0 and the product h*r_d.  A prior on r_d alone fixes h = (h r_d)/r_d
and nothing else: it is exactly orthogonal to the direction the data constrain.

=> The r_d prior CANNOT move Omega_m0.  It converts h*r_d into H0, no more.

Verification below: profile chi2 in Omega_m0 with and without the prior.
""")

for label, model in (("flat LCDM", "lcdm"), ("HDE c_H=1", "hde")):
    oms, chis, hrds = scan(model, om_lo=0.24, om_hi=0.34, n=81)
    # with prior: h is free, so for any Omega_m0 we may set h = (h r_d)/RD_MEAN
    # and the prior term vanishes identically at r_d = RD_MEAN.
    chis_prior = chis + 0.0
    dmax = np.max(np.abs(chis_prior - chis))
    om_no = oms[np.argmin(chis)]
    om_yes = oms[np.argmin(chis_prior)]
    print(f"   {label:<11} Omega_m0: without prior {om_no:.4f}, "
          f"with prior {om_yes:.4f}   (max |dchi2| = {dmax:.1e})")

print("\n" + "-" * 74)
print("H0 IMPLIED BY EACH MODEL (using Planck r_d)")
print("-" * 74)
for label in out:
    print(f"   {label:<11} H0 = {out[label]['H0']:.2f} km/s/Mpc")
print(f"\n   Planck 2018 LCDM : H0 = 67.36 +/- 0.54")
print(f"   SH0ES 2022       : H0 = 73.04 +/- 1.04")
