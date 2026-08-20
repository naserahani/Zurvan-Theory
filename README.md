# Zurvan Theory

**Naser Ahani**
Department of Physics, Amirkabir University of Technology, Tehran, Iran
naser.ahani@aut.ac.ir

Official repository for the Zurvan research programme: sources, figures, and numerical code for papers on emergent spacetime, holographic information, and the foundations of physics.

---

## The framework in brief

Zurvan proposes that spacetime, matter, and forces are not primary entities but **emergent patterns** of a single underlying oscillatory entity. Two commitments follow:

- **Time as rhythm, not axis.** Time is not a pre-existing dimension along which events are arranged; it is the name we give to the succession of fundamental events. Irreversibility is then structural rather than something to be explained away.
- **Information as the bridge to observation.** The framework's contact with data runs through horizon thermodynamics and the holographic principle, which is where its quantitative predictions come from.

The name is taken from the Iranian mythological figure of boundless time. The framework is a research programme in progress, not a finished theory; open problems are listed below and stated explicitly in each paper.

---

## Papers

### Dark Energy as Inverse Holographic Information (2026) — *flagship*
`papers/informational-lambda/`

The cosmological constant problem is reformulated in informational terms via the exact identity ρ_Λ = (3/8) ρ_p / N_H, where N_H is the Bekenstein–Hawking information content of the horizon. Promoting this identity from a statement about the asymptotic horizon to a dynamical law referred to the instantaneous future event horizon yields holographic dark energy with the usually free parameter **fixed to c = 1** by the one-bit-per-4ℓ_p² normalization.

The result is a prediction with **no adjustable parameters**:

> **w₀ = −1/3 − (2/3)√Ω_Λ0 ≈ −0.885**

Evaluated with the Planck Ω_Λ0 = 0.685 this gives w₀ ≈ −0.885. A direct fit of the model to the DESI DR2 BAO distance ratios reproduces the data as well as flat ΛCDM (Δχ² = 0.07 with the same number of parameters); see `calculations/p6-desi-fit/`. The model also predicts a definite thawing history whose sign of evolution differs from the current CPL central fit, which makes it sharply falsifiable by forthcoming data. The paper states clearly what the framework does *not* explain: the value of Λ itself remains an input, and the explanatory burden is relocated rather than dissolved.

### A Hypothesis for a Sub-Planckian Timescale (2025)
`papers/holographic-timescale/`

Derives the informational constant N_Z ≈ 3.3 × 10¹²² from the holographic entropy of the cosmic horizon and the associated timescale T_pp = T_p / N_Z. Interpretive rather than derivational: the number is the known de Sitter entropy; what the framework adds is a reading of it as a count of fundamental events per Planck time.

### A Scalar–Tensor Model for Singularity Resolution (2025)
`papers/singularity-scalar-tensor/`

A scalar field non-minimally coupled to the matter Lagrangian undergoes a density-triggered phase transition, driving the total stress–energy to zero and replacing the classical singularity with a locally flat core, while reducing exactly to general relativity at low density. **Work in progress** — see `STATUS.md` in that directory for the open technical gaps; not ready for journal submission.

### On the Reversibility of Time (2026)
`papers/time-reversal-essay/`

A philosophy-of-physics essay arguing that debates over time-reversal operators presuppose time as an external parameter, and examining what changes if that presupposition is dropped.

---

## Open problems

The framework's unresolved questions, in priority order:

1. **Bell's theorem.** Any claim of sub-quantum determinism must confront Bell. The route is nonlocal hidden variables, in the tradition of Bohm and 't Hooft — locality, not determinism, is what the theorem rules out.
2. **Lorentz invariance.** A sequential process defines a preferred frame. Either macroscopic invariance must be recovered as emergent, or a suppressed violation consistent with observational bounds must be predicted.
3. **The memory problem.** If the fundamental entity is at one location per fundamental step, what carries the state of everywhere else between visits?
4. **The phase rule.** Deriving the Feynman weight from the underlying dynamics, rather than restating it. Until this is done, the framework's reading of the path integral remains an internal document rather than a publishable result.
5. **Completing the singularity model.** See `papers/singularity-scalar-tensor/STATUS.md`.
6. **Direct data fit.** Fitting the parameter-free distance–redshift relation to DESI BAO and supernova compilations.

## What this framework no longer claims

Earlier presentations included claims now withdrawn as incompatible with established physics or as overstated:

- A distance-based mechanism for the fundamental forces — incompatible with QCD confinement and with the gauge structure of the Standard Model.
- The elimination of quantum uncertainty as such — pending item 1 above.
- A *solution* to the cosmological constant problem — the informational reading reformulates it.

Superseded material is kept in `retired/` rather than deleted.

---

## Repository layout

```
papers/         sources, figures, and compiled PDFs, one directory per paper
calculations/   numerical scripts supporting the papers
retired/        superseded versions and withdrawn claims
```

## Citing

Please cite the individual paper. Preprints are also available on ResearchGate and Zenodo.

## Contact

Correspondence, criticism, and collaboration are welcome: naser.ahani@aut.ac.ir
