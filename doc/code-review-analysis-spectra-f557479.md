# Code review: analysis workflow and `singleParticleSpectra` at f557479

## Scope

This review covers `src/main.cpp`, `src/Analysis.{h,cpp}`, and `src/single_particleSpectra.{h,cpp}` at f557479. It starts from the state after the fixes in 3244b93, which are described in `doc/analysis-spectra-investigation.md`, and does not repeat those findings. Each candidate below was checked either with synthetic events of known content or by tracing the code path.

## Findings

| # | Severity | Finding |
|---|---|---|
| 1 | High | φ(1020) yields in the flow workflow are 2.04× too large |
| 2 | Medium | Error bars for three- and four-particle correlators ignore fluctuations in the tuple count |
| 3 | Low, pre-existing | The ensemble script expects rapidity-correlator filenames that the workflow does not write |
| 4 | Low | `check_dNdSV` writes out of bounds when a smearing width exceeds the bin width |

### 1. φ(1020) yields are divided by a branching ratio that does not apply

For `particle_monval == 333`, `singleParticleSpectra` divides the following by `reconst_branching_ratio = 0.489`:

- the integrated yield and its error ([output_Qn_vectors](../src/single_particleSpectra.cpp#L760))
- the pT spectrum ([differential spectrum](../src/single_particleSpectra.cpp#L831))
- dN/dy ([rapidity density](../src/single_particleSpectra.cpp#L2301))

This rescaling assumes every φ was reconstructed from a K⁺K⁻ pair. In the flow workflow that never happens:

- `run_flow_spectra` builds one `particleSamples` from the global configuration ([reader](../src/Analysis.cpp#L127)). That reader sets `reconst_flag` only when the global `particle_monval` is 333 ([constructor](../src/particleSamples.cpp#L118)), and the default is 211.
- Modes 0, 9, and 10 clear `reconst_flag` when reading ([reader](../src/particleSamples.cpp#L404)).

The φ entries in the flow output are therefore the explicitly stored φ records, divided by 0.489.

**Evidence.** The input was 200 events in `read_in_mode=9`, each with one stored φ at y ~ U(−0.4, 0.4) and 20 pions, run with default flow settings. The expected dN/dy in |y| < 0.5 is 1.0.

| Output | Reported | Expected |
|---|---|---|
| `particle_333_vndata_y_-0.5_0.5.dat`, row 0 | 2.04498978 | 1.0 |
| `particle_333_dNdy_pT_0.2_3.dat`, bins inside \|y\| < 0.4 | ≈ 2.6 per unit y | ≈ 1.25 per unit y |

**Suggested fix.** Apply the branching ratio only to φ produced by `perform_particle_reconstruction`. For example, the reader could report whether reconstruction ran. Explicit φ records should keep weight 1.

When reconstruction does run (global `particle_monval=333` with a UrQMD, JAM, or SMASH reader), stored φ records are counted alongside the reconstructed ones. A mixed sample would then need separate weights.

### 2. Multiparticle error bars ignore fluctuations in the tuple count

The correlators reported are R = ΣX_e/ΣY_e, where X_e is the event's tuple numerator and Y_e its ordered tuple count. Their errors are computed as √(Var(X)/N)/⟨Y⟩ in each of these writers:

- [Cmnk](../src/single_particleSpectra.cpp#L1702)
- [Cmnk rapidity](../src/single_particleSpectra.cpp#L1886)
- [Cn4](../src/single_particleSpectra.cpp#L2089)
- [SCmn](../src/single_particleSpectra.cpp#L2127)
- the same-sign and opposite-sign variants

This formula leaves out the covariance between X and Y. In 3244b93, `vn{2}` was moved to the delta method, √Σ(X_e − R Y_e)²/ΣY_e, which accounts for it. The higher-order correlators still use the old formula.

The term the formula misses grows like R·σ(Y), which scales as M⁴ for four-particle correlators, while the sampling noise grows like M². The error is therefore overestimated in wide centrality bins at high multiplicity. Accurate error bars for a non-zero four-particle correlator need Σ X·Y and Σ Y² accumulated alongside the existing sums.

**Evidence.** Independent replicas of synthetic events were generated with a random event plane, v₂ = 0.1, |η| < 1, and pT ~ U(0.25, 1.5) GeV. The table compares the reported error, averaged over replicas, with the standard deviation of the reported values across replicas. The delta-method error was computed in Python from the same events.

| Multiplicity M per event | Correlator | Reported err / replica scatter | Delta method / replica scatter |
|---|---|---|---|
| U{100..400}, 150 events, 60 replicas | C₂{4} numerator ⟨4⟩ | 1.05 | 1.06 |
| U{100..400}, 150 events, 60 replicas | C₁₁₂ ⟨3⟩ | 0.98 | 0.92 |
| U{1000..3000}, 100 events, 40 replicas | C₂{4} numerator ⟨4⟩ | **1.85** | 1.14 |
| U{1000..3000}, 100 events, 40 replicas | C₁₁₂ ⟨3⟩ (≈ 0) | 0.99 | 0.96 |

With 40 replicas, the scatter itself is uncertain by about 11%. In the high-multiplicity case, the central values were unchanged: ⟨4⟩ = 1.01×10⁻⁴ against v₂⁴ = 1.0×10⁻⁴.

These errors propagate downstream: `analysis_multi_particle_correlation_within_UrQMDevent_h5.py` averages the SCmn error column ([line 443](../ebe_scripts/analysis_multi_particle_correlation_within_UrQMDevent_h5.py#L443)).

**Suggested fix.** Accumulate Σ X·Y and Σ Y² for every correlator, as `vn{2}` already does, and reuse `write_pair_moments`.

### 3. Ensemble script filenames do not match the workflow (pre-existing)

`analysis_multi_particle_correlation_within_UrQMDevent_h5.py` reads `particle_9999_Cmnk_%s_pT_0.2_3.dat` ([line 304](../ebe_scripts/analysis_multi_particle_correlation_within_UrQMDevent_h5.py#L304)). The correlation definitions in `FlowAnalysis` and `FlowAnalysis_LHC` set the reference pT to 0.2–2.0 GeV ([Analysis.cpp](../src/Analysis.cpp#L348), [LHC](../src/Analysis.cpp#L826)), so the workflow writes `..._pT_0.2_2.dat`. The verification run produced `particle_9999_Cmnk_eta12_pT_0.2_2.dat`.

The 0.2–2.0 GeV setting was already present at 6a63746, so this predates the recent commits. Either the script or the definitions should be updated, and the script should fail clearly when a file is missing.

### 4. `check_dNdSV` can write past its arrays

The τ, x, y, and η_s histograms compute `idx = (value + random − min)/width` with `random ∈ [0, intrinsic_*)` and no upper-bound check ([τ](../src/single_particleSpectra.cpp#L2385), [x](../src/single_particleSpectra.cpp#L2397), [η_s](../src/single_particleSpectra.cpp#L2458)). Only the r histogram checks its index.

The last bin overflows when the smearing exceeds one bin width. With the hard-coded grids that means:

- `intrinsic_dtau` > 0.294 fm
- `intrinsic_dx` > 0.408 fm
- `intrinsic_detas` > 0.154

The defaults (0.01, 0.1, 0.1) are safe. The problem only arises with `check_spatial_dis=1`. This was derived from the code and not executed.

**Suggested fix.** Add the same `idx < N` guard the r histogram uses.

## Checked and not reported

- **`_weakFD` labels in modes 0, 9, and 10.** The reader clears its `resonance_weak_feed_down_flag` member there, but the decayer chooses `pdg_weak.dat` when it is constructed ([particle_decay.cpp](../src/particle_decay.cpp#L51)). The cleared member has no other reader. Weak decays still run, so the labels are correct; the reset is dead code.
- **The second, non-weak-feed-down flow pass reruns every definition.** The old code at 6a63746 also reran the full flow analysis, so this is unchanged behavior. `FlowAnalysis_multistrange_particles` is not called from the workflow.
- **Three-particle overlap subtractions.** For the integrated charge-dependent flags 0, 1, 3, and 4, the subtraction terms were derived by hand and match the code. The kernels were enumeration-tested in 3244b93.
- **Rapidity-difference bin offset.** `static_cast<int>(rapidity_dis_min / drap)` is correct only when that ratio is an integer. The workflow hard-codes the −2…2, 41-point grid, for which −2/0.1 is exactly −20.
- **Rapidity abscissa.** It uses (grid value + Σy)/(N + 1), deliberately seeded with the grid value. This only matters in sparse bins.

## Reproduction

- **Build and run.** `build/src/hadronic_afterburner_tools.e read_in_mode=9 analyze_flow=1 randomSeed=12345`, using the root `parameters.dat` and `EOS/`.
- **Event format.** Each event is a native `int` count, followed by `struct.pack("=i9f", pdg, m, t, x, y, z, E, px, py, pz)` per track.
- **Generator.** Python with numpy, `default_rng(seed)`:
  - M ~ U{Mmin..Mmax} pions alternating +211 and −211, m = 0.13957 GeV
  - pT ~ U(0.25, 1.5) GeV, η ~ U(−1, 1)
  - φ sampled by rejection from 1 + 2v₂ cos 2(φ − Ψ), with a random Ψ per event
  - optionally, φ(1020) records with PDG 333, m = 1.019 GeV, pT = 0.8 GeV, y ~ U(−0.4, 0.4)
- **Finding 1 input.** 200 events with M = 20, v₂ = 0, and one φ each (seed 1).
- **Finding 2 inputs.** Replica i used seed 100 + i:
  - 60 replicas of 150 events, M ∈ [100, 400], v₂ = 0.08
  - 40 replicas of 100 events, M ∈ [1000, 3000], v₂ = 0.1
- **Finding 2 extraction.** Row 2 of `particle_9999_Cn4_eta_-1_1.dat` and row 1 of `particle_9999_Cmnk_eta_-1_1.dat`. In Python, the per-event X_e for n = 2 is |Q₂|⁴ + |Q₄|² − 2Re(Q₄Q₂*²) − 4(M−2)|Q₂|² + 2M(M−3), with Y_e = M(M−1)(M−2)(M−3). For C₁₁₂, X_e is Re(Q₁Q₁Q₂*) − |Q₂|² − 2|Q₁|² + 2M, with Y_e = M(M−1)(M−2).

No source files were changed during this review.

**Status.** Findings 1 and 2 are fixed; see `doc/fix-phi-branching-ratio-and-correlator-errors.md`.
