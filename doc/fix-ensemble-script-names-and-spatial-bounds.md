# Fix: ensemble correlation script filenames and spatial histogram bounds

This fixes findings 3 and 4 of `doc/code-review-analysis-spectra-f557479.md`, starting from ac01e17.

## 3. `analysis_multi_particle_correlation_within_UrQMDevent_h5.py` filenames

### Problem

The mismatch was wider than the review reported. With `compute_correlation = 1`, `FlowAnalysis` and `FlowAnalysis_LHC` write the charged-hadron correlators ([Analysis.cpp](../src/Analysis.cpp#L333)) as follows:

- **Pseudorapidity windows:** |η| < 1 and |η| < 2, giving `particle_9999_*_eta_-1_1.dat` and `*_eta_-2_2.dat`.
- **Reference particles:** 0.2 < pT < 2 GeV, so the rapidity-difference files are `*_pT_0.2_2.dat`.

The script read two sets of names that the workflow never writes:

- `_eta_-0.5_0.5` for `vn2`, `Cn2_ss/os`, `Cmnk(_ss/_os)`, `Cn4`, `SCmn`, and `vndata_diff`
- `_pT_0.2_3` for `vn2_eta12(_ss/_os)` and `Cmnk_eta12/eta13(_ss/_os)`

`h5py`'s `.get()` returns `None` for a missing dataset. On a real ensemble, the script therefore stopped at its first correlator with `IndexError: too many indices for array`, which does not name the missing file.

The script was added in 2018–2019. The C++ workflow that writes these files was added in March 2020 (90f56f4) and has used these windows ever since. The script is the stale side, so the workflow definitions are unchanged.

### Change

- `rap_region = "-1_1"`, which matches the |η| < 1 window and the script's STAR output names. `ref_pT_region = "0.2_2"` is used for the rapidity-difference files. Both sit at the top of the script with a comment pointing to the workflow definitions.
- Every per-event read goes through `load()`. It exits with status 1 and names the HDF5 file, the event group, the missing dataset, and the flags that produce it.
- Estimators, column indexing, and output files are unchanged.

## 4. Out-of-bounds writes in `check_dNdSV`

### Problem

The τ, x, y, η_s, and (τ, x) histograms used `idx = (value + smear − min)/width` with `smear ∈ [0, intrinsic_*)` and no upper bound. Bin widths are (max − min)/(N − 1), so a smear wider than one bin can index past the last bin. The threshold widths are 0.294 fm for τ, 0.408 fm for x and y, and 0.154 for η_s. The defaults (0.01, 0.1, 0.1) are safe. This only affects `check_spatial_dis = 1`.

### Change

- Every histogram now checks `0 ≤ idx < N` before writing, as the r histogram already did. Entries smeared past the last bin are dropped, also as the r histogram already did.
- The check comes after the random draw, so the random-number sequence and the output for in-range smearing are unchanged.
- Negative indices could not occur (value > min, smear ≥ 0); the lower check is defensive.

## Validation

### Ensemble script

The fixture was four "hydro events", each a separate run of 50 events from the scratch generator described in the review document. The settings were:

- M ∈ U{200..400} pions, v₂ = 0.08, random Ψ, |η| < 1, pT ~ U(0.25, 1.5) GeV
- generator seeds 201–204
- `read_in_mode=9 analyze_flow=1 randomSeed=12345` with the root `parameters.dat`, which has `compute_correlation`, `flag_charge_dependence`, and `compute_corr_rap_dep` set to 1

The `particle_9999_*` outputs were packed with `ebe_scripts/combine_results_into_hdf5.py` (h5py 3.16, numpy 2.5).

| Check | Before | After |
|---|---|---|
| Run on the ensemble | `IndexError` at the first correlator | Exit 0, all 15 output files written |
| Layouts of the files read (`vn2` 10×5, `vn2_eta12` 41×41, `Cmnk` 9×3, `Cmnk_ss` 9×5, `Cmnk_eta12` 41×19, `Cn4` 5×3, `SCmn` 6×3) | n/a | Match the script's column indices |
| C₁₁₂ (charged) | n/a | 1.51374×10⁻⁶; direct count-weighted mean over events: 1.51374×10⁻⁶ |
| v₂{2}² | n/a | 5.55821×10⁻³; direct: 5.55821×10⁻³ |
| C₂{4} | n/a | −3.87169×10⁻⁵; direct: −3.87169×10⁻⁵. v₂{4} = 0.0789 against input 0.08 |
| One dataset removed from one event | n/a | Exit 1: "… has no dataset particle_9999_Cmnk_eta13_os_pT_0.2_2.dat in event /event_3 …" |

### Spatial histograms

The new unit test `Analysis.SpatialHistogramsDropSmearsPastTheLastBin` uses 400 tracks just inside the upper τ, x, and y − η_s edges, with all smearing widths set to 1.

| Build | Result |
|---|---|
| ac01e17 with AddressSanitizer | `heap-buffer-overflow` in `check_dNdSV` |
| Fixed, with AddressSanitizer | Passes, no report |

The kept fractions match the analytic acceptance for a uniform smear:

| Histogram | Kept | Analytic acceptance |
|---|---|---|
| τ | 150/400 | 39% |
| x | 204/400 | 51% |
| η_s | 72/400 | 20% |

All 35 unit tests pass (20 particleSamples, 4 HBT, 11 Analysis). The build has no warnings, and `git clang-format` was applied to the changed C++ lines.

## Limitations and related issues not fixed

- **Sibling scripts.** Two other scripts read the same stale names. I did not change them:
  - `analysis_multi_particle_correlation_within_UrQMDevent_h5_minimumbias.py` already uses `-1_1`, but reads `_pT_0.2_3` for the rapidity-difference files.
  - `analysis_multi_particle_correlation_within_UrQMDevent.py` (text input) reads `-0.5_0.5` and `_pT_0.2_3`.
- **Symmetric-cumulant error (pre-existing).** In the SC section of the fixed script, `vm2_ch_err` subtracts `vn2_ch**2`, a stale and already normalized value from the previous loop iteration, where `vm2_ch**2` is meant. This affects only the error column of `symmetric_cumulant_STAR.dat`.
- **Script estimators.** The script's ensemble estimators combine per-event errors as `Σ w·err/Σ w/√nev`. They were not reviewed here.
