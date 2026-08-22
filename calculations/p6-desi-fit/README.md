# P6 — The $c_H = 1$ holographic dark energy model against DESI DR2 BAO and Pantheon+

Numerical support for *Dark Energy as Inverse Holographic Information*
(`papers/informational-lambda/`).

## The model

Holographic dark energy with the **future event horizon** as the infrared cutoff:

$$\rho_\Lambda = \frac{3 c_H^2 M_p^2}{R_h^2}, \qquad R_h = a\!\int_t^\infty \frac{dt'}{a(t')}$$

The parameter $c_H$, normally free and fitted, is **fixed to 1** by the
Bekenstein–Hawking one-bit-per-$4\ell_p^2$ normalisation. This leaves

$$\frac{d\Omega_\Lambda}{d\ln a} = \Omega_\Lambda(1-\Omega_\Lambda)\left(1 + \frac{2\sqrt{\Omega_\Lambda}}{c_H}\right), \qquad w_\Lambda(z) = -\frac13 - \frac23\frac{\sqrt{\Omega_\Lambda(z)}}{c_H}$$

**Free parameters: $\Omega_{m0}$ and $h\,r_d$ — two, the same count as flat ΛCDM.**

## Data

**DESI DR2 Results II, Table IV** — 13 independent BAO measurements over
$0.295 \le z \le 2.330$: isotropic $D_V/r_d$ for BGS, and $D_M/r_d$ with $D_H/r_d$ for six
anisotropic tracers. $\chi^2$ uses a 2×2 covariance block per anisotropic tracer, built
from the published $\sigma_M$, $\sigma_H$ and correlation $\rho_{M,H}$. Correlations
*between* redshift bins are not modelled. Degrees of freedom: $13-2 = 11$.

**Pantheon+** (Brout et al. 2022) — 1701 supernovae in the release, **1580 used** after
$z_{\rm HD} > 0.01$ and removal of calibrators. STAT+SYS covariance. The absolute
magnitude (and hence $H_0$) is absorbed analytically into a single nuisance offset, so the
SN likelihood constrains the *shape* of the distance–redshift relation only.

### Getting the Pantheon+ files

They are ~33 MB and are **not committed** to this repository. Download them into a
`pantheon/` subdirectory next to the scripts:

```bash
mkdir -p pantheon && cd pantheon
B="https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/Pantheon%2B_Data/4_DISTANCES_AND_COVAR"
curl -L "$B/Pantheon%2BSH0ES.dat"          -o data.txt
curl -L "$B/Pantheon%2BSH0ES_STAT%2BSYS.cov" -o cov.txt
```

## Pipeline validation

Before fitting anything, the same code was run on flat ΛCDM and compared with published
values:

| | this code | published |
|---|---|---|
| $\Omega_{m0}$ (BAO) | 0.2973 | 0.2975 ± 0.0086 (DESI DR2) |
| $h\,r_d$ [Mpc] | 101.55 | 101.54 ± 0.73 (DESI DR2) |
| $\Omega_{m0}$ (SN) | 0.3326 | 0.334 ± 0.018 (Brout et al. 2022) |

## Result 1 — BAO

| model | $\Omega_{m0}$ | $h\,r_d$ [Mpc] | $\chi^2$ | $\chi^2/\mathrm{dof}$ |
|---|---|---|---|---|
| flat ΛCDM | 0.2973 ± 0.0085 | 101.55 ± 0.72 | 10.539 | 0.958 |
| HDE $c_H=1$ | 0.2717 ± 0.0084 | 99.72 ± 0.73 | 10.605 | 0.964 |

$\Delta\chi^2 = +0.066$ with the same number of parameters. Releasing $c_H$ gives
$c_H = 0.925\,(-0.129\,{+}0.163)$: **the fixed value $c_H=1$ lies 0.46σ from the free best
fit**, and fixing it costs $\Delta\chi^2 = 0.24$.

## Result 2 — supernovae, and SN/BAO concordance

Fitting each probe separately within each model:

| model | SN only | BAO only | combined | SN vs BAO |
|---|---|---|---|---|
| flat ΛCDM | 0.3326 (−0.0180 +0.0185) | 0.2973 | 0.3042 | **1.75σ** |
| HDE $c_H=1$ | 0.2722 (−0.0157 +0.0163) | 0.2717 | 0.2718 | **0.03σ** |

**The two independent probes agree within the $c_H=1$ model and disagree within ΛCDM.**
This is the substantive finding: the model's preference for a lower $\Omega_{m0}$ is not a
BAO artefact — the supernovae, which carry no sound horizon and no ΛCDM-derived input,
independently prefer the same value.

On the combined data, $\Delta\chi^2 = -3.45$ in favour of the model, at equal parameter
count.

## Result 3 — the $\omega_m$ tension

BAO fixes $h\,r_d$; the model fixes $\Omega_{m0}$; and $r_d$ itself depends on
$\omega_m = \Omega_m h^2$. Solving the three self-consistently (`omh2.py`):

| | $\omega_m$ | $h$ | $r_d$ [Mpc] |
|---|---|---|---|
| flat ΛCDM | 0.1405 | 0.6875 | 147.69 |
| HDE $c_H=1$ | **0.1114** | 0.6402 | 155.80 |
| Planck 2018 | **0.1430 ± 0.0011** | — | — |

**This is the model's outstanding difficulty, and it is stated here rather than left for a
reader to find.** $\omega_m$ is fixed by the acoustic peak structure and damping tail at
$z \approx 1100$ — early-universe physics that late-time dark energy does not alter. In
this model $\Omega_\Lambda$ at recombination is $\simeq 1.7\times10^{-4}$, changing $H$ by
under 0.01% there, so the CMB constraint applies to it as well.

*No σ value is quoted for this difference.* The Planck number is inferred within ΛCDM;
quoting a σ between it and a non-ΛCDM model's requirement would be the same
methodological error as comparing a single parameter across model families. The
statement is qualitative and deliberately so: **the model requires an $\omega_m$ well
below the CMB value.** Settling it requires a joint analysis with a CMB likelihood, not a
number computed here.

## A structural note on $r_d$ priors

Every DESI BAO observable scales as $1/(h\,r_d)$. Writing $u = c/(100\,h\,r_d)$, the model
vector is $u\cdot m(\Omega_{m0})$ and $\chi^2$ is exactly quadratic in $u$. **A prior on
$r_d$ is therefore exactly orthogonal to the direction BAO constrains**: it converts
$h\,r_d$ into $H_0$ and moves $\Omega_{m0}$ by nothing at all. `p6_1_fast.py` verifies this
numerically — the profile $\chi^2$ in $\Omega_{m0}$ is bit-identical with and without the
prior.

This is why the supernova test, not an $r_d$ prior, was the informative one.

## Running

```bash
pip install numpy scipy matplotlib
python3 p6_fit.py         # BAO fits; writes p6_results.npy
python3 p6_errors.py      # 1-sigma intervals, c_H constraint, residuals, figure
python3 p6_1_fast.py      # closed-form BAO scan; the r_d prior result
python3 omh2.py           # self-consistent omega_m (no data needed)
python3 wz_history.py     # Table 1 and Figure 1 of the paper
python3 p6_2_pantheon.py  # supernovae + combination   [needs pantheon/]
```

Run from this directory: `p6_errors.py`, `p6_1_fast.py` and `p6_2_pantheon.py` all import
from `p6_fit.py`. Outputs are written next to the scripts.

`p6_2_pantheon.py` takes a few seconds. It uses continuous minimisation
(`minimize_scalar`, bounded Brent, `xatol=1e-5`) rather than a grid: the true minimum
never lands on a grid node, and a 0.005-step grid is enough to shift the ΛCDM SN/BAO
tension from 1.75σ to 1.96σ. The `OMS` grid that remains in the file is used only to draw
the figure.

Tested with numpy 1.26, scipy 1.11, matplotlib 3.6.

## Files

| file | contents |
|---|---|
| `p6_fit.py` | model ODE, distance observables, $\chi^2$, BAO fits |
| `p6_errors.py` | 1σ intervals, $c_H$ constraint, per-tracer residuals, figure |
| `p6_1_fast.py` | closed-form BAO scan; demonstrates the $r_d$-prior result |
| `p6_2_pantheon.py` | Pantheon+ likelihood, SN/BAO concordance, combination |
| `omh2.py` | self-consistent $\omega_m$, $h$, $r_d$ per model |
| `wz_history.py` | integrates the model from the fitted $\Omega_{\Lambda 0}$; writes `wz_table.tex` and `wz_figure.pdf`, which are **Table 1 and Figure 1 of the paper** |
| `wz_table.tex` | the eight table rows, as generated |
| `p6_hde_vs_desi.png` | $w(z)$ history and BAO distance ratios |
| `p6_2_sn_bao.png` | SN, BAO and combined $\Omega_{m0}$ likelihoods per model |
| `wz_figure.pdf` | $w(z)$ and $\Omega_\Lambda(z)$ histories — Figure 1 of the paper |

## Caveats

- **No CMB likelihood.** That is exactly where the $\omega_m$ question is decided.
- **Profile likelihood, not MCMC.** Adequate for two parameters; a publication fit should
  use MCMC.
- **$r_d$ from the DESI fitting formula**, evaluated at $\omega_m \approx 0.111$, which is
  22% from where that formula was calibrated. The extrapolation error is ~1%.
- **Radiation neglected** in $E(z)$; ~0.05% at $z = 2.33$ against ~1.2% data precision.
