# P6 — Direct fit of the $c_H = 1$ holographic dark energy model to DESI DR2 BAO

Numerical support for *Dark Energy as Inverse Holographic Information*
(`papers/informational-lambda/`).

## The model

Holographic dark energy with the **future event horizon** as the infrared cutoff:

$$\rho_\Lambda = \frac{3 c_H^2 M_p^2}{R_h^2}, \qquad R_h = a\!\int_t^\infty \frac{dt'}{a(t')}$$

The parameter $c_H$, normally free and fitted, is **fixed to 1** by the
Bekenstein–Hawking one-bit-per-$4\ell_p^2$ normalisation. This leaves

$$\frac{d\Omega_\Lambda}{d\ln a} = \Omega_\Lambda(1-\Omega_\Lambda)\left(1 + \frac{2\sqrt{\Omega_\Lambda}}{c_H}\right), \qquad w_\Lambda(z) = -\frac13 - \frac23\frac{\sqrt{\Omega_\Lambda(z)}}{c_H}$$

**Free parameters: $\Omega_{m0}$ and $h\,r_d$ — two, the same count as flat ΛCDM.**

## The data

**DESI DR2 Results II, Table IV** — 13 independent measurements over
$0.295 \le z \le 2.330$:

- BGS at $z=0.295$: isotropic $D_V/r_d$
- Six anisotropic tracers (LRG1, LRG2, LRG3+ELG1, ELG2, QSO, Lyα): $D_M/r_d$ and
  $D_H/r_d$ per tracer

$\chi^2$ uses a **2×2 covariance block per anisotropic tracer**, built from the published
$\sigma_M$, $\sigma_H$ and the correlation coefficient $\rho_{M,H}$, plus the scalar BGS
term. Correlations *between* redshift bins are not modelled — the published covariance is
not available in that form. Degrees of freedom: $13 - 2 = 11$.

## Pipeline validation

Before fitting the model, the same code was run on flat ΛCDM and compared with the
official DESI numbers:

| | this code | DESI DR2 official |
|---|---|---|
| $\Omega_{m0}$ | 0.2973 | 0.2975 ± 0.0086 |
| $h\,r_d$ [Mpc] | 101.55 | 101.54 ± 0.73 |

Agreement to four significant figures. The pipeline reproduces the published result.

## Main result

| model | $\Omega_{m0}$ | $h\,r_d$ [Mpc] | $\chi^2$ | $\chi^2/\mathrm{dof}$ |
|---|---|---|---|---|
| flat ΛCDM | 0.2973 ± 0.0085 | 101.55 ± 0.72 | 10.539 | 0.958 |
| HDE $c_H=1$ | 0.2717 ± 0.0084 | 99.72 ± 0.73 | 10.605 | 0.964 |

$$\Delta\chi^2 = +0.066 \quad \text{with the same number of parameters}$$

The model describes the DESI BAO distance ratios as well as flat ΛCDM does.

**Testing $c_H = 1$.** Releasing $c_H$ gives $c_H = 0.925\,(-0.129\,{+}0.163)$ with
$\chi^2 = 10.365$. **The fixed value $c_H = 1$ lies 0.46σ from the free best fit**, and
fixing it costs $\Delta\chi^2 = 0.24$. The data exert no pressure on this postulate.

Equation of state at the best fit: $w_0 = -0.902$, rising toward $-1/3$ with redshift
(thawing, effective $w_a = +0.311$ matched at $z=0.5$).

## Running

```bash
pip install numpy scipy matplotlib
python3 p6_fit.py      # fits LCDM, HDE(c_H=1), HDE(c_H free); writes p6_results.npy
python3 p6_errors.py   # uncertainties, c_H constraint, residuals; writes the figure
```

`p6_errors.py` imports from `p6_fit.py`, so run them from this directory in that order.
Both write their output next to the scripts. Total runtime is a few minutes; the cost is
the future-horizon integral, which is evaluated out to $z_{\max}$ for every likelihood
call.

Tested with numpy 1.26, scipy 1.11, matplotlib 3.6.

## Files

| file | contents |
|---|---|
| `p6_fit.py` | model ODE, distance observables, $\chi^2$, profile-likelihood fits |
| `p6_errors.py` | 1σ intervals, $c_H$ constraint, per-tracer residuals, figure |
| `p6_hde_vs_desi.png` | $w(z)$ history and BAO distance ratios against the data |

## Caveats

- **BAO only.** No CMB likelihood, no supernovae. The model prefers
  $\Omega_{m0} = 0.2717 \pm 0.0084$, which is lower than the ΛCDM value inferred from the
  same BAO data and lower still than the Planck value; adding CMB would test that
  directly.
- **Profile likelihood, not MCMC.** Adequate for two parameters; a publication fit should
  use MCMC.
- **Radiation neglected** in $E(z)$. At $z = 2.33$ the effect is ~0.05%, against data
  precision of ~1.2%.
