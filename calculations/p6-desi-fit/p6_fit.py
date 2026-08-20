#!/usr/bin/env python3
"""
P6 — Direct fit of the c_H = 1 holographic dark energy model (Zurvan M2)
against DESI DR2 BAO, compared with flat LCDM.

Model:
    rho_Lambda = 3 c_H^2 M_p^2 / R_h^2,  R_h = future event horizon,  c_H = 1 (FIXED)
    => dOmega_L/dlna = Omega_L (1-Omega_L)(1 + 2 sqrt(Omega_L)/c_H)
    => w_L(z)        = -1/3 - (2/3) sqrt(Omega_L(z))/c_H

Free parameters: Omega_m0 and (h * r_d).  Same count as flat LCDM.
"""

import numpy as np
from pathlib import Path
from scipy.integrate import solve_ivp, quad
from scipy.interpolate import interp1d
from scipy.optimize import minimize

C_KMS = 299792.458
C_OVER_H100 = C_KMS / 100.0          # = 2997.92458 Mpc  (so c/H0 = C_OVER_H100/h)

# ---------------------------------------------------------------- DESI DR2 data
# Values: DESI DR2 Results II, Table IV.  13 independent measurements.
BGS = dict(z=0.295, DV=7.942, sV=0.075)

ANISO = [
    dict(name="LRG1",      z=0.510, DM=13.588, sM=0.167, DH=21.863, sH=0.425, rho=-0.459),
    dict(name="LRG2",      z=0.706, DM=17.351, sM=0.177, DH=19.455, sH=0.330, rho=-0.404),
    dict(name="LRG3+ELG1", z=0.934, DM=21.576, sM=0.152, DH=17.641, sH=0.193, rho=-0.416),
    dict(name="ELG2",      z=1.321, DM=27.601, sM=0.318, DH=14.176, sH=0.221, rho=-0.434),
    dict(name="QSO",       z=1.484, DM=30.512, sM=0.760, DH=12.817, sH=0.516, rho=-0.500),
    dict(name="Lya",       z=2.330, DM=38.988, sM=0.531, DH= 8.632, sH=0.101, rho=-0.431),
]

for t in ANISO:
    cov = np.array([[t["sM"]**2, t["rho"]*t["sM"]*t["sH"]],
                    [t["rho"]*t["sM"]*t["sH"], t["sH"]**2]])
    t["icov"] = np.linalg.inv(cov)

Z_MAX = 2.6

# ---------------------------------------------------------------- HDE solver
def hde_E(Om0, cH=1.0, zmax=Z_MAX, npts=4000):
    """Return E(z)=H/H0 interpolator and Omega_Lambda(z) interpolator for HDE."""
    OL0 = 1.0 - Om0

    def rhs(x, y):
        OL = np.clip(y[0], 1e-12, 1 - 1e-12)
        return [OL * (1 - OL) * (1 + 2 * np.sqrt(OL) / cH)]

    x_end = np.log(1.0 / (1.0 + zmax))          # x = ln a, integrate backwards
    xs = np.linspace(0.0, x_end, npts)
    sol = solve_ivp(rhs, (0.0, x_end), [OL0], t_eval=xs,
                    rtol=1e-10, atol=1e-12, method="RK45")
    if not sol.success:
        raise RuntimeError("ODE failed")

    x = sol.t
    OL = sol.y[0]
    a = np.exp(x)
    z = 1.0 / a - 1.0

    # flat, matter + HDE:  1 - OL = Om0 a^-3 / E^2
    E = np.sqrt(Om0 * a**-3 / (1.0 - OL))

    idx = np.argsort(z)
    return (interp1d(z[idx], E[idx], kind="cubic", bounds_error=False,
                     fill_value="extrapolate"),
            interp1d(z[idx], OL[idx], kind="cubic", bounds_error=False,
                     fill_value="extrapolate"))


def lcdm_E(Om0):
    return lambda z: np.sqrt(Om0 * (1 + z)**3 + (1 - Om0))


# ---------------------------------------------------------------- observables
def observables(Efun, hrd, zlist):
    """Return dict z -> (DM/rd, DH/rd, DV/rd)."""
    K = C_OVER_H100 / hrd                    # = (c/H0)/rd
    zs = np.linspace(0.0, Z_MAX, 3000)
    integ = np.concatenate([[0.0], np.cumsum(np.diff(zs) * 0.5 *
                            (1.0 / Efun(zs[:-1]) + 1.0 / Efun(zs[1:])))])
    chi = interp1d(zs, integ, kind="cubic")

    out = {}
    for z in zlist:
        DM = K * float(chi(z))
        DH = K / float(Efun(z))
        DV = (DM**2 * z * DH) ** (1.0 / 3.0)
        out[z] = (DM, DH, DV)
    return out


def chi2(Efun, hrd):
    zlist = [BGS["z"]] + [t["z"] for t in ANISO]
    obs = observables(Efun, hrd, zlist)

    c2 = ((obs[BGS["z"]][2] - BGS["DV"]) / BGS["sV"]) ** 2
    for t in ANISO:
        d = np.array([obs[t["z"]][0] - t["DM"], obs[t["z"]][1] - t["DH"]])
        c2 += float(d @ t["icov"] @ d)
    return c2


# ---------------------------------------------------------------- fitting
def fit_hde(cH=1.0):
    def nll(p):
        Om0, hrd = p
        if not (0.05 < Om0 < 0.7 and 60 < hrd < 160):
            return 1e6
        try:
            E, _ = hde_E(Om0, cH)
        except Exception:
            return 1e6
        return chi2(E, hrd)
    r = minimize(nll, [0.30, 101.0], method="Nelder-Mead",
                 options=dict(xatol=1e-6, fatol=1e-8, maxiter=4000))
    return r


def fit_lcdm():
    def nll(p):
        Om0, hrd = p
        if not (0.05 < Om0 < 0.7 and 60 < hrd < 160):
            return 1e6
        return chi2(lcdm_E(Om0), hrd)
    r = minimize(nll, [0.30, 101.0], method="Nelder-Mead",
                 options=dict(xatol=1e-6, fatol=1e-8, maxiter=4000))
    return r


if __name__ == "__main__":
    N_DATA = 13
    N_PAR = 2
    DOF = N_DATA - N_PAR

    print("=" * 68)
    print("P6 — c_H = 1 holographic dark energy vs DESI DR2 BAO")
    print("=" * 68)

    rl = fit_lcdm()
    Om_l, hrd_l = rl.x
    print(f"\nflat LCDM        chi2 = {rl.fun:8.3f}  (dof={DOF})  "
          f"chi2/dof = {rl.fun/DOF:5.3f}")
    print(f"                 Om0  = {Om_l:.4f}   h*rd = {hrd_l:.2f} Mpc")

    rh = fit_hde(1.0)
    Om_h, hrd_h = rh.x
    OL0 = 1 - Om_h
    w0 = -1.0/3.0 - (2.0/3.0) * np.sqrt(OL0)
    print(f"\nHDE  c_H = 1     chi2 = {rh.fun:8.3f}  (dof={DOF})  "
          f"chi2/dof = {rh.fun/DOF:5.3f}")
    print(f"                 Om0  = {Om_h:.4f}   h*rd = {hrd_h:.2f} Mpc")
    print(f"                 w0   = {w0:.4f}  (from Omega_L0 = {OL0:.4f})")

    print(f"\nDelta chi2 (HDE - LCDM) = {rh.fun - rl.fun:+.3f}   "
          f"(same number of parameters)")

    # free-c_H reference
    def nll_c(p):
        Om0, hrd, cH = p
        if not (0.05 < Om0 < 0.7 and 60 < hrd < 160 and 0.3 < cH < 3.0):
            return 1e6
        try:
            E, _ = hde_E(Om0, cH)
        except Exception:
            return 1e6
        return chi2(E, hrd)
    rc = minimize(nll_c, [0.30, 101.0, 1.0], method="Nelder-Mead",
                  options=dict(xatol=1e-6, fatol=1e-8, maxiter=6000))
    print(f"\nHDE  c_H free    chi2 = {rc.fun:8.3f}  (dof={N_DATA-3})")
    print(f"                 Om0  = {rc.x[0]:.4f}   h*rd = {rc.x[1]:.2f}   "
          f"c_H = {rc.x[2]:.4f}")

    # w(z) history for the c_H=1 best fit
    _, OLf = hde_E(Om_h, 1.0)
    print("\n  w(z) history, c_H = 1 best fit:")
    for z in [0.0, 0.3, 0.5, 1.0, 1.5, 2.0, 2.33]:
        OL = float(OLf(z))
        w = -1.0/3.0 - (2.0/3.0)*np.sqrt(OL)
        print(f"    z = {z:4.2f}   Omega_L = {OL:.4f}   w = {w:+.4f}")

    # CPL-equivalent w_a from a two-point match at z=0 and the pivot
    OL_p = float(OLf(0.5))
    w_p = -1.0/3.0 - (2.0/3.0)*np.sqrt(OL_p)
    wa = (w_p - w0) / (0.5/1.5)
    print(f"\n  effective w_a (matched at z=0.5): {wa:+.4f}")

    np.save(Path(__file__).with_name("p6_results.npy"),
            dict(lcdm=(Om_l, hrd_l, rl.fun),
                 hde1=(Om_h, hrd_h, rh.fun, w0),
                 hdec=(rc.x[0], rc.x[1], rc.x[2], rc.fun),
                 dof=DOF), allow_pickle=True)
