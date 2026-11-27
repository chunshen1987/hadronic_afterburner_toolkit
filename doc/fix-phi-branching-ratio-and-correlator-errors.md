# Fix: φ(1020) branching ratio and multiparticle correlator errors

This fixes findings 1 and 2 of `doc/code-review-analysis-spectra-f557479.md`, starting from f557479.

## 1. Branching ratio applied only to reconstructed φ

### Problem

For `particle_monval == 333`, every φ yield was divided by 0.489, the K⁺K⁻ branching ratio. In the flow workflow no φ is ever reconstructed, so the stored φ records came out 2.045× too high. When reconstruction did run (UrQMD, JAM, or SMASH reader with global `particle_monval=333`), stored φ were rescaled along with the rebuilt ones.

### Change

- `particle_info` gains `bool reconstructed = false`. `particleSamples::perform_particle_reconstruction` sets it on each φ it builds from a K⁺K⁻ pair. No reader, writer, or `sizeof` depends on the struct layout.
- `singleParticleSpectra` counts reconstructed tracks separately. A reconstructed track has weight w = 1/0.489 and a stored track has weight 1:
  - **Integrated yield** (`*_vndata_*`, row 0 and row 99): N = N_all + (w − 1)·N_rec. The Poisson error is √(N_all + (w² − 1)·N_rec)/(N_ev·Δy), from the sum of squared weights.
  - **pT spectrum** (`*_vndata_diff_*`): the same rule in each pT bin, using the spectrum's rapidity window and pT bins.
  - **dN/dy** (`*_dNdy_*`): the same rule in each rapidity bin, with no pT cut, as for the existing counts.
  - **Spatial checks** (`check_333_*`, only with `check_spatial_dis=1`): one scale factor, 1 + (w − 1)·N_rec/N, measured on the tracks that enter the histograms. This is exact when the sample is all reconstructed or all stored, and approximate for a mixture.
- Mean pT and vₙ are unweighted, as before.
- Particles other than φ have w = 1 and are unaffected.

## 2. Delta-method errors for three- and four-particle correlators

### Problem

The correlators are R = ΣX_e/ΣY_e, where X_e is the event's tuple numerator and Y_e its ordered tuple count, which is row or entry 0. Their errors were √(Var X/N)/⟨Y⟩, which ignores the X–Y covariance. In wide multiplicity bins this overestimated the C₂{4} numerator error by 1.85×.

### Change

- For each correlator, the code now also accumulates Σ X_e·Y_e and reports the delta-method error of a ratio estimator under independent events:

  σ_R = √(ΣX² − 2R·ΣXY + R²·ΣY²) / ΣY = √Σ(X_e − R·Y_e)² / ΣY

  This is the estimator `vn{2}` has used since 3244b93.
- Bins with ΣY ≤ 0 report 0 for both value and error, as before.
- Row 0 (the mean tuple count and its standard error) is unchanged.
- It applies to:
  - `Cmnk` and the charge-dependent `Cmnk_ss` and `Cmnk_os`, both columns 12 and 13
  - `Cn4` and `SCmn`
  - `Cmnk_eta12` and `Cmnk_eta13`, and their `_ss` and `_os` variants
- Estimators, reference particles (`vn_rapidity_dis_pT_min/max`, the `rap_min/max` window), and self-correlation subtraction are unchanged. Event loops now pass zeroed per-event buffers to the existing kernels, and a shared helper accumulates Σ, Σ², and the cross term.
- Kernel signatures are unchanged. The two charge-dependent drivers take optional cross-moment pointers.
- File names, headers, and column layouts are unchanged. The four writers share helpers (`tuple_count`, `tuple_ratio`, `write_four_particle_correlation`).

## Validation

**Unit tests.** All 34 pass (20 particleSamples, 4 HBT, 10 Analysis), run from a copy of `unit_tests/` with `EOS/`. There are two new regression tests in `unit_tests/Analysis_unittest.cc`. Both fail on f557479 and pass with the fix.

| Test | Setup | f557479 | Fixed |
|---|---|---|---|
| `ReconstructedPhiAloneCarriesTheBranchingRatio` | One UrQMD event with a stored φ and a K⁺K⁻ pair from a φ decay; checks integrated dN/dy, its error, row 99, the pT-spectrum sum, and the dN/dy integral | 4.090 (both φ ÷ 0.489) | 3.04499 = 1 + 1/0.489 |
| `MultiparticleErrorsUseTheRatioDeltaMethod` | Three mode-9 events of 4, 5, and 7 pions; X_e and Y_e enumerated over ordered distinct tuples for C₂{4} ⟨4⟩ and C₁₁₂ | error 0.0478 | 0.0368, matches enumeration |

**End-to-end runs** (`read_in_mode` as listed, root `parameters.dat`, `analyze_flow=1`, `randomSeed=12345`):

| Input | Quantity | f557479 | Fixed | Expected |
|---|---|---|---|---|
| Mode 9: 200 events, one stored φ each, \|y\| < 0.4 | φ dN/dy, \|y\| < 0.5 | 2.04499 | 1.00000 | 1.0 |
| Mode 1: 300 events, each with one φ → K⁺K⁻ at a shared vertex and one stored φ, \|y\| < 0.3 | φ dN/dy, \|y\| < 0.5 | 4.08998 | 3.04499 | 1 + 1/0.489 |

**Output comparison.** Fixed and f557479 outputs were byte-compared on three inputs: the two φ runs above and 100 events with M ∈ U{1000..3000}, v₂ = 0.1 (136–137 files each).
- The only files that differ are the φ yield files and the multiparticle correlator files.
- In the correlator files, headers, value columns, and count rows are identical. Only the error columns change.

**Replicas.** There were 40 independent replicas of 100 events each, with M ∈ U{1000..3000}, v₂ = 0.1, random Ψ, |η| < 1, and seeds 101–140.

| Correlator | Mean | Replica scatter | Reported error, f557479 | Reported error, fixed |
|---|---|---|---|---|
| C₂{4} ⟨4⟩ | 1.008×10⁻⁴ | 6.59×10⁻⁶ | 1.85× scatter | 7.59×10⁻⁶ (1.15×) |
| C₁₁₂ | −1.9×10⁻⁷ | 4.78×10⁻⁶ | 0.99× | 4.73×10⁻⁶ (0.99×) |

- On replicas 1–20, the reported errors equal a Python delta-method calculation from the same events to 6 significant digits.
- With 40 replicas, the scatter itself is uncertain by about 11%.

**Build and format.** The build has no warnings, and `git clang-format` was applied to the changed lines.

## Limitations

- The spatial check histograms use one sample-wide weight, not a per-bin one, as noted above.
- The delta method assumes independent events. Oversamples of the same collision should be grouped before being treated as independent, as for `vn{2}`.
- The event generators used for the end-to-end and replica checks were scratch scripts and are not in the repository. The event format and generator settings are described in the review document.
