# Verification of code-review findings 1–4 on revision 3244b93

## Question

A code review of the `dev` branch at 3244b93 reported four candidate problems:

1. HBT internal mixing uses only the preceding collision, which may be a non-independent oversample.
2. The HBT pair loops no longer preselect tracks, which makes the O(N²) loop more expensive.
3. Balance-function internal mixing can pair an event with itself, which adds self-pairs at Δy = 0.
4. Σ0 → Λγ feed-down draws a new decay on every `filter_particles` call.

Each claim was tested with synthetic events whose answer is known. Confirmed problems were fixed.

## Summary

| # | Verdict | Change |
|---|---------|--------|
| 1 | Not confirmed: the radii and λ are unbiased in every mixing mode tested | None |
| 2 | Confirmed as a runtime cost. The acceptance change was intentional | Lossless track preselection; 2.5–6× faster on wide-rapidity input, byte-identical output |
| 3 | Confirmed: like-sign mixed Δy = 0 bin inflated about 2.9× with 5-event batches | Same-batch mixing draws only from other events |
| 4 | Confirmed: all 12 Σ0 daughters differed between lists of one batch | One cached decay per event record and track |

## 1. HBT mixing reference

**Claim.** The default mixed reference is the preceding collision group, after a random azimuthal rotation. The first collision gets no mixed pairs. When oversamples of one hydro event are read as separate collisions (`number_of_oversample_events = 1`), the reference comes from the same source.

**Findings.** The first collision does enter the denominator: when collision 2 is processed, its pairs with collision 1 are binned. Only the first and last collisions have one neighbour instead of two. The mixed-pair count was about 1.9 times the same-event count in every internal-mixing run.

The bias test used a static Gaussian source with R = 4 fm along each axis, t = 0, and positions independent of momenta. For that source, the cosine-weighted estimator gives C − 1 = exp(−q²R²/ħc²) with λ = 1. Pion momenta carried v₂ = 0.15 about a plane shared by the oversamples of one hydro event. Fits of λ exp(−q²R²/ħc²) to |q| shells below 0.1 GeV, using the toolkit normalization (`pair_ratio` × mixed), gave:

| Input and mixing | λ | R (fm) | same/mixed at \|q\| > 60 MeV |
|---|---|---|---|
| 300 independent events, internal (previous collision) | 1.022 ± 0.014 | 4.025 ± 0.029 | 0.9982 |
| same events, external stream of 300 independent events | 1.010 ± 0.014 | 4.013 ± 0.029 | 0.9987 |
| 30 hydro × 10 oversamples, group size 1 (oversample as reference) | 0.995 ± 0.014 | 3.966 ± 0.029 | 0.9996 |
| 30 hydro × 10 oversamples, group size 10 | 0.995 ± 0.004 | 3.985 ± 0.009 | 0.9997 |

All four runs recover the input radius. λ is within 1.6σ of 1. Rotating the reference removes the shared reaction plane, so an oversample of the same hydro event is a valid reference for this estimator. The earlier code also mixed within one file, so it used the same kind of reference. No code change was made.

**Limitation.** This source has no position–momentum correlation. In a flowing source, the difference between mixing within one hydro event and across hydro events is second order in event-by-event fluctuations. It was not quantified here. When every record forms a single collision group, internal mixing has no reference, as README.md already notes.

## 2. HBT track preselection

**Claim.** Commit 3169027 removed the single-track `tanh(HBTrap)` preselection, so every track enters the pair loops.

**Findings.** Applying `HBTrap_min/max` to the pair rapidity was intentional. It was finding 3 of `doc/code-review-9ce7e81.md` and is documented in README.md, so that change stays. The cost is real: on wide-rapidity input, every track enters both pair loops.

**Fix.** `src/HBT_correlation.cpp` now removes, before pairing, only tracks that cannot belong to any binned pair:

- Every binned pair has q_inv² < 3Q², with Q = max(|q_min − Δq/2|, |q_max + Δq/2|). The 3D bins need |q_out|, |q_side|, |q_long| < Q, and q_inv² ≤ q_out² + q_side² + q_long² in the lab or the LCMS. The q_inv bins need q_inv < Q.
- q_inv² = 2p₁·p₂ − m₁² − m₂² ≥ 2m_T1 m_T2 (cosh Δy − 1) − (m₁ − m₂)². This uses p_T1·p_T2 ≤ p_T1 p_T2 and m_T1 m_T2 − p_T1 p_T2 ≥ m₁m₂ (Cauchy–Schwarz). Together with the first point, this bounds |Δy| < Δy_max for each track, using the smallest partner m_T.
- K_z/K_E is an energy-weighted mean of the two track velocities, so the pair rapidity lies between y₁ and y₂. A track is kept only if HBTrap_min − Δy_max < y < HBTrap_max + Δy_max.
- A 10⁻³ relative margin absorbs rounding. If any track is not timelike, or has E ≤ |p_z|, every track in that call is kept.

Kept tracks stay in their original order, so the pair order and the `needed_number_of_pairs` caps are unchanged. The pair kinematics, the estimator, and the normalization are unchanged.

**Evidence.** HEAD was built in a separate worktree for the baseline. Times are wall-clock times for the whole run:

| Input | HEAD | Fixed | Output files identical |
|---|---|---|---|
| 3 collisions × 20000 π, \|y\| < 4, group size 10 | 15.5 s | 2.48 s | 2/2 |
| 20 events × 2000 π, \|y\| < 4, n_KT = 5, qnpts = 31 | 1.64 s | 0.66 s | 8/8 |

In the heavy case, each collision kept about 6400 of 20000 tracks. Byte-for-byte `cmp` of all `HBT_*` outputs also passed for these configurations:

- `invariant_radius_flag = 1`: 12 files
- `azimuthal_flag = 1`: 16 files
- each run in the table for issue 1, including external mixing

## 3. Balance-function self-mixing

**Claim.** With `read_in_real_mixed_events = 0`, the mixed lists are filtered from the same `full_particle_list`. The random draw of `iev_mixed` can therefore equal `iev`.

**Evidence.** The input had 2000 events, each with one π⁺ and one π⁻ at independent y ~ U(−1, 1) and p_T = 0.5 GeV. It was run with `event_buffer_size = 10`, which gives 5-event batches. In the mixed histogram, the share of like-sign pairs in the central Δy bin (|Δy| < 0.1) was:

| Build | like-sign (SS) | unlike-sign (OS) | expected without self-mixing |
|---|---|---|---|
| HEAD | 0.2807 | 0.0925 | 0.0975 |
| Fixed | 0.0978 | 0.0943 | 0.0975 |

Self-mixing with probability 1/5 predicts 0.8 × 0.0975 + 0.2 = 0.278 for SS, which matches HEAD.

**Fix.** `combine_and_bin_mixed_particle_pairs` takes a `same_batch` flag. `calculate_balance_function` sets it when no external mixed stream is used. In that case the reference event is drawn uniformly from the other nev − 1 events of the batch. One random draw per event is kept, as before. A batch with a single event has no independent reference, so its mixed pairs are skipped and a warning is printed once. Previously those pairs were self-mixed. Pair counting, the Δy and Δφ binning, the p_T cuts, and the normalization `(same/N_same)/(mixed/N_mixed)` are unchanged.

## 4. Σ0 → Λγ feed-down consistency

**Claim.** `filter_particles` decays each Σ0 again on every call.

**Evidence.** The input had 4 events with three Σ0 each. It was run with `particle_alpha = 3122` and `particle_beta = −3122`, which makes the `bbar` list a Λ list too, plus `particle_monval = 3122` and the Σ flag on. On HEAD, the BF a list differed from the bbar list, from the a_mixed list, and from the bbar_mixed list for all 12 tracks. After the fix, all 12 match.

**Fix.** `particleSamples::get_sigma_daughter` caches the Λ daughter per event record and track index. Every list filtered from the same batch therefore sees the same decay: the PoI list, the anti list, the BF lists, and the same-batch mixed lists. Cache entries are erased in `clear_out_previous_record`, which is the only place event records are freed. The Λ or anti-Λ sign follows the Σ0 sign. The feature still applies only when a filter selects ±3122, as the parameter comment says.

The review's side remark that the feature "is active for any configured `particle_monval`" describes intended behavior: a BF Λ list also receives Σ0 feed-down.

## Checks performed

- New regression tests:
  - `BalanceFunction.SameBatchMixingNeverPairsAnEventWithItself`, which covers 20 repeats and the single-event batch.
  - `particleSamples.SigmaFeedDownSharesOneDecayAcrossLists`.
- Both new tests fail on HEAD and pass with the fix. All 32 unit tests pass (20 + 4 + 8). They were run from a copy of `unit_tests/` with `EOS/` beside it.
- An AddressSanitizer build (Debug, `CMAKE_CXX_FLAGS_DEBUG`) produced no reports in these runs:
  - a 300-event Σ0 run with 7-track buffers (43 batches), with internal and external mixing
  - HBT with group sizes 1 and 10, the q_inv branch, and the azimuthal branch
- No compiler warnings with `-Wall -Wextra`. The changed lines were formatted with `git clang-format --style=file HEAD` (clang-format 21.1.8). Untouched lines in the edited files were left as they were. After formatting, the tests and the HBT byte comparison against HEAD passed again.

## Unresolved limitations

- The preselection bound uses the smallest partner m_T for every track, so it is conservative. Sorting by rapidity would prune more, but it would change the pair order and therefore which pairs are kept under `needed_number_of_pairs`.
- Pre-existing review item 5 was not changed. When the external mixed stream ends before the main file, the run aborts and writes no output. One ASan run hit this when the mixed file was shorter than the main file.
- For issue 1, no comparison with a flowing source that has position–momentum correlations was made.

## Reproduction

The scripts below live in the session scratchpad and are reproduced here. All runs use `read_in_mode = 9`, the root `parameters.dat` and `EOS/`, and `randomSeed = 12345` unless noted.

Binary record per event: a native `int` count, then `struct.pack("=i9f", pdg, m, t, x, y, z, E, px, py, pz)` per track.

- **HBT inputs:** `gen_hbt.py out nhydro noversample npart seed [R=4] [v2=0.15] [ymax=0.6]`. p_T ~ |N(0.35, 0.12)| GeV, y ~ U(−ymax, ymax), φ sampled from 1 + 2v₂ cos 2(φ − Ψ), and x, y, z ~ N(0, R) fm.
  - Issue 1 used `indep.bin` (300 1 300 seed 1), `indep_mixed.bin` (seed 2), and `overs.bin` (30 10 300 seed 3).
  - Issue 2 used `wide.bin` (20 1 2000 seed 11, ymax 4) and `heavy.bin` (3 10 2000 seed 21, ymax 4).
- **HBT command:** `analyze_flow=0 analyze_HBT=1 long_comoving_boost=0 qnpts=21 q_min=-0.1 q_max=0.1 n_KT=2 KT_min=0.15 KT_max=0.55 HBTrap_min=-0.5 HBTrap_max=0.5`.
  - Add `read_in_real_mixed_events=1` or `number_of_oversample_events=10` as listed.
  - The timing runs override `long_comoving_boost=1`; the wide case also overrides `qnpts=31 q_min=-0.15 q_max=0.15 n_KT=5`.
- **HBT fit:** read `HBT_raw_correlation_function_KT_*.dat`. Recover `pair_ratio` as the processed-file denominator divided by raw mixed counts. Form C − 1 = Σcos / (pair_ratio × mixed) in 10 shells of q² up to 0.01 GeV², then fit λ exp(−q²R²/ħc²) with ħc = 0.19733 GeV fm.
- **Balance function:** 2000 events, each with one π⁺ and one π⁻, y ~ U(−1, 1), p_T = 0.5 GeV, numpy seed 7. Run `analyze_balance_function=1 particle_alpha=211 particle_beta=-211 event_buffer_size=10`. Compare column 8 (SS mixed) and column 5 (OS mixed) of `Balance_function_211_-211_Delta_y.dat` in the central row.
- **Σ0:** 4 events with three Σ0 each (PDG 3212, m = 1.192642 GeV, p_T = 0.8 GeV, y ~ U(−1, 1), numpy seed 8). A C++11 driver linked to `libhadronic_afterburner_tools` calls `read_in_particle_samples_and_filter()` and `read_in_particle_samples_mixed_event_and_filter()`. It then compares BF lists by `event_track_id` with `particle_monval=3122 particle_alpha=3122 particle_beta=-3122 resonance_weak_feed_down_Sigma_to_Lambda_flag=1`.
