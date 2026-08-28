#!/usr/bin/env python3
"""
P6.3 -- Joint analysis: CMB distance priors + DESI DR2 BAO + Pantheon+.

This is the analysis the reviewers asked for, in its compressed form.

Why distance priors rather than the full Planck likelihood.  The reviewers'
objection was that correlated shifts in h, omega_b, n_s, A_s and tau might
absorb the omega_m discrepancy.  The Planck distance priors (R, l_A, omega_b)
are obtained from the Planck 2018 chains *after marginalising over exactly
those parameters*, so using them addresses that objection directly.  Chen,
Huang & Wang (2019) verified against LCDM, wCDM and CPL that the priors
reproduce the full-likelihood constraints.  The present model modifies only the
late-time background (Omega_Lambda ~ 2e-4 at recombination), which is the
regime where the compression is valid.

What this is NOT: a full Boltzmann + likelihood analysis.  The compression
assumes standard pre-recombination physics and inherits a residual dependence
on the LCDM chains from which the priors were derived.

Parameters: Omega_m0, h, omega_b.  SN absolute magnitude marginalised
analytically; r_d computed from the same sound-horizon integral as r_s(z_*).

Data:
  CMB   Chen, Huang & Wang 2019 (arXiv:1808.05724), Planck TT,TE,EE+lowE
  BAO   DESI DR2 Results II, 13 measurements
  SN    Pantheon+ (Brout et al. 2022), STAT+SYS
"""

from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.interpolate import interp1d
from scipy.optimize import minimize

from p6_fit import BGS, ANISO, C_OVER_H100

# numpy renamed trapz -> trapezoid in 2.0; support both.
_trapz = getattr(np, "trapezoid", None) or np.trapz

C_KMS = 299792.458
TCMB = 2.7255
NEFF = 3.046
FT = (TCMB / 2.7) ** -4                      # (T/2.7K)^-4

# ---------------------------------------------------------------- CMB priors
CMB_D = np.array([1.750235, 301.4707, 0.02235976])
CMB_ICOV = np.array([[9.43923971e4, -1.3604913e3, 1.6645172916e6],
                     [-1.3604913e3,  1.614349e2,  3.6716180e3],
                     [1.6645172916e6, 3.6716180e3, 7.97191825162e7]])

# ---------------------------------------------------------------- background
def omega_r(om_b_unused=None):
    """Physical radiation density (photons + massless neutrinos)."""
    og = 2.4728e-5 * (TCMB / 2.7255) ** 4
    return og * (1.0 + 0.2271 * NEFF)

OM_R = omega_r()


def E_lcdm_factory(Om0, h):
    orh2 = OM_R
    Or0 = orh2 / h**2
    OL = 1.0 - Om0 - Or0
    return lambda z: np.sqrt(Om0 * (1 + z)**3 + Or0 * (1 + z)**4 + OL)


def E_hde_factory(Om0, h, cH=1.0, zmax=3000.0, npts=24000):
    """HDE with radiation.  Exact ODE:
         dOmega_L/dx = Omega_L [ (1-Omega_L)(1 + 2 sqrt(Omega_L)/cH) + Omega_r ]
       with Omega_r = Or0 a^-4 (1-Omega_L) / (Om0 a^-3 + Or0 a^-4).

    Above zmax the dark component is utterly negligible (Omega_L ~ a, so
    ~1e-7 by z=1e6), and E reverts exactly to the matter+radiation form.
    The returned callable uses the ODE solution below zmax and that limit
    above it, so the sound-horizon integral (which reaches a ~ 1e-8) is
    evaluated correctly rather than by extrapolation.
    """
    Or0 = OM_R / h**2
    OL0 = 1.0 - Om0 - Or0

    def rhs(x, y):
        OL = float(np.clip(y[0], 1e-30, 1 - 1e-12))
        a = np.exp(x)
        Or = Or0 * a**-4 * (1 - OL) / (Om0 * a**-3 + Or0 * a**-4)
        return [OL * ((1 - OL) * (1 + 2 * np.sqrt(OL) / cH) + Or)]

    x_end = np.log(1.0 / (1.0 + zmax))
    xs = np.linspace(0.0, x_end, npts)
    s = solve_ivp(rhs, (0.0, x_end), [OL0], t_eval=xs, rtol=1e-10, atol=1e-16)
    if not s.success:
        raise RuntimeError("HDE ODE failed")
    a = np.exp(s.t)
    z = 1.0 / a - 1.0
    E = np.sqrt((Om0 * a**-3 + Or0 * a**-4) / (1.0 - s.y[0]))
    i = np.argsort(z)
    lo = interp1d(z[i], E[i], kind="cubic", bounds_error=False, fill_value=np.nan)

    def Efun(zq):
        zq = np.asarray(zq, dtype=float)
        hi = np.sqrt(Om0 * (1 + zq)**3 + Or0 * (1 + zq)**4)
        out = np.where(zq <= zmax, lo(np.clip(zq, 0.0, zmax)), hi)
        return out if out.ndim else float(out)

    return Efun


# ---------------------------------------------------------------- sound horizon
def z_star(om_m, om_b):
    """Hu & Sugiyama 1996 photon-decoupling redshift."""
    g1 = 0.0783 * om_b**-0.238 / (1 + 39.5 * om_b**0.763)
    g2 = 0.560 / (1 + 21.1 * om_b**1.81)
    return 1048 * (1 + 0.00124 * om_b**-0.738) * (1 + g1 * om_m**g2)


def z_drag(om_m, om_b):
    """Eisenstein & Hu 1998 drag epoch."""
    b1 = 0.313 * om_m**-0.419 * (1 + 0.607 * om_m**0.674)
    b2 = 0.238 * om_m**0.223
    return 1291 * om_m**0.251 / (1 + 0.659 * om_m**0.828) * (1 + b1 * om_b**b2)


def r_s(zend, Efun, h, om_b):
    """Comoving sound horizon, c/H0 * int_0^a da/(a^2 E sqrt(3(1+R_b)))."""
    a_end = 1.0 / (1.0 + zend)
    Rb_coeff = 31500.0 * om_b * FT

    def integ(a):
        z = 1.0 / a - 1.0
        return 1.0 / (a**2 * float(Efun(z)) * np.sqrt(3.0 * (1 + Rb_coeff * a)))

    val, _ = quad(integ, 1e-8, a_end, limit=200, epsabs=1e-10, epsrel=1e-9)
    return (C_OVER_H100 / h) * val


def chi_comoving(z, Efun, h):
    zs = np.linspace(0.0, z, 4000)
    Ev = Efun(zs)
    return (C_OVER_H100 / h) * _trapz(1.0 / Ev, zs)


# ---------------------------------------------------------------- likelihoods
def chi2_cmb(Efun, Om0, h, om_b):
    om_m = Om0 * h**2
    zst = z_star(om_m, om_b)
    DM = chi_comoving(zst, Efun, h)
    R = np.sqrt(Om0) * h * 100.0 * DM / C_KMS
    lA = np.pi * DM / r_s(zst, Efun, h, om_b)
    d = np.array([R, lA, om_b]) - CMB_D
    return float(d @ CMB_ICOV @ d), R, lA


def chi2_bao(Efun, h, rd):
    K = 1.0 / rd
    zs = np.linspace(0.0, 2.6, 3000)
    integ = np.concatenate([[0.0], np.cumsum(np.diff(zs) * 0.5 *
                            (1 / Efun(zs[:-1]) + 1 / Efun(zs[1:])))])
    chi = interp1d(zs, integ, kind="cubic")
    pre = C_OVER_H100 / h

    def DM(z): return pre * float(chi(z)) * K
    def DH(z): return pre / float(Efun(z)) * K

    z = BGS["z"]
    dv = (DM(z)**2 * z * DH(z)) ** (1 / 3)
    out = ((dv - BGS["DV"]) / BGS["sV"])**2
    for t in ANISO:
        d = np.array([DM(t["z"]) - t["DM"], DH(t["z"]) - t["DH"]])
        out += float(d @ t["icov"] @ d)
    return out


if __name__ == "__main__":
    print("=" * 72)
    print("  P6.3 -- CMB distance priors + DESI DR2 BAO + Pantheon+")
    print("=" * 72)
    print(f"\nradiation: omega_r = {OM_R:.6e}")
    # sanity: LCDM Planck best fit
    Om0, h, ob = 0.3153, 0.6736, 0.02237
    E = E_lcdm_factory(Om0, h)
    c2, R, lA = chi2_cmb(E, Om0, h, ob)
    rd = r_s(z_drag(Om0 * h**2, ob), E, h, ob)
    print(f"\nLCDM at Planck best fit (Om={Om0}, h={h}, ob={ob}):")
    print(f"   z_*  = {z_star(Om0*h**2, ob):.2f}   z_d = {z_drag(Om0*h**2, ob):.2f}")
    print(f"   R    = {R:.4f}   (Planck {CMB_D[0]:.4f})")
    print(f"   l_A  = {lA:.3f}   (Planck {CMB_D[1]:.3f})")
    print(f"   r_d  = {rd:.2f} Mpc   (Planck 147.09)")
    print(f"   chi2_CMB = {c2:.2f}")
