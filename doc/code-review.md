# Scientific code review

Reviewed revision `6a63746` for defects that can bias observables, interrupt analysis, or produce invalid output. The review sampled the event readers, flow and correlation estimators, resonance decays, and Python ensemble processing. Line references in the findings refer to the reviewed revision.

## High priority

### 1. The three-particle correlation versus Δη₁₃ uses the wrong harmonic

The `flag == 2` branch integrates particle 2 at harmonic `k = n + m` and uses that sum in the raw triplet term. The estimator requires `Q_m` for particle 2, and the overlap correction for particles 1 and 3 also uses that integrated sum. This affects the inclusive and charge-dependent Δη₁₃ paths. In a synthetic event with three distinct particles at identical azimuth `π/3` and in the same rapidity bin, the six ordered distinct triplets each contribute `cos(φ + φ − 2φ) = 1`. The code returned `C000 = 6` and `C112 = −3` in that bin, giving `−0.5` after triplet normalization instead of `1`. Use harmonic `m` for the integrated vector and validate the overlap terms against explicit enumeration for inclusive, same-sign, and opposite-sign samples. [Harmonic accumulation and use](../src/single_particleSpectra.cpp#L1657) · [Inclusive caller](../src/single_particleSpectra.cpp#L415) · [Charge-dependent callers](../src/single_particleSpectra.cpp#L1574)

### 2. SMASH binary header parsing can exceed buffer bounds

The header reader uses a file-provided string length for a fixed-size buffer without checking capacity. It also passes the four bytes of the magic field to C-string operations without initializing a terminator. Even a valid magic field can therefore be read past its end. Check each header read, compare the magic as four bytes, and read the version into bounded storage with explicit length handling. The same reader also shares one format-version member between the primary and mixed streams, so different input versions can make the primary reader skip the wrong number of bytes at each event footer. Keep version metadata per stream or reject mismatched streams explicitly. This finding is from source inspection. [Stream setup](../src/particleSamples.cpp#L155) · [Header reader](../src/particleSamples.cpp#L288) · [Version-dependent event footer](../src/particleSamples.cpp#L1116)

### 3. HBT HDF5 averaging cannot consume current C++ output

The azimuthally integrated C++ output has five columns with `ecoOutput=0` and two with `ecoOutput=1`. The HDF5 averager expects a legacy table with at least eight, including raw same-event and mixed-event pair counts. It reads column 5 as the mixed-event count, while the current five-column output ends at column 4 and can already have a pair normalization applied. Moving columns alone would lose information needed for ensemble normalization. Define a common format that retains the raw pair counts and weighted numerator, then validate the producer and consumer together. A scratch HDF5 file containing current executable output reproduced an out-of-range column read after a test-only change in memory let the script get past centrality selection. [C++ output writer](../src/HBT_correlation.cpp#L738) · [HDF5 accumulator](../ebe_scripts/average_event_HBT_correlation_function_h5.py#L174)

### 4. NumPy shadows the scalar `min` call in several centrality scripts

The scripts import `*` from NumPy and call `min(a, b)` with two integer bounds. In the installed NumPy version this resolves to `numpy.min(a, axis=b)`. On a 20-event input the HBT averaging script failed with an `AxisError` before reading its correlation tables. The same import and call pattern appears in the minimum-bias flow, photon, spin, and multiparticle scripts. Use the built-in operation explicitly or switch to qualified NumPy imports. Also check that a percentile scan includes the least-multiplicitous event in the last class and has a defined policy for tied multiplicities. The existing `> lower_cut` test excludes events tied to the minimum when that value is used as the 100% edge. [HBT](../ebe_scripts/average_event_HBT_correlation_function_h5.py#L134) · [Flow](../ebe_scripts/average_event_spvn_h5_minimumbias.py#L1597) · [Photon](../ebe_scripts/average_event_photon_h5_minimumbias.py#L305) · [Spin](../ebe_scripts/average_event_spin_h5_minimumbias.py#L1065) · [Multiparticle](../ebe_scripts/analysis_multi_particle_correlation_within_UrQMDevent_h5_minimumbias.py#L87)

## Medium priority

### 5. Resonances at rest produce nonfinite daughter momenta

Both two- and three-body decays use `(γ − 1) (v · p) / |v|²` in the boost, which evaluates to `0/0` for a mother at rest. The two-body routine reproduced this for `ρ⁰ → π⁺π⁻`, returning `nan` in all three daughter momentum components despite a finite energy of `0.3845 GeV`. The boost coefficient can be evaluated as `γ²/(γ+1)` times `v · p`, which has a finite rest limit. Cover an at-rest and a moving mother with a four-momentum conservation check. [Two-body boost](../src/particle_decay.cpp#L305) · [Three-body boost](../src/particle_decay.cpp#L438)

### 6. Rapidity shifts change the HBT phase by transforming only momentum

`boostParticles` shifts `E` and `pz` without shifting `t` and `z`, but HBT evaluates `q⁰Δt − qᶻΔz` using both. For `q⁰ = 0`, `qᶻ = 0.2 GeV`, `Δt = 0`, `Δz = 1 fm`, and no transverse contribution, a common longitudinal boost by `0.7` leaves the cosine at approximately `0.529`. Boosting only momentum changes it to `0.294` for the same pair. This is separate from an intended change in acceptance after the shift. Transform both quantities consistently for HBT, or evaluate the phase from the original four-vectors and specify the frame of the cuts. The numerical example checks the phase algebra independently of pair selection. [Momentum shift](../src/particleSamples.cpp#L442) · [Phase calculation](../src/HBT_correlation.cpp#L431)

### 7. HBT radii fitting can fail at grid reshape and at output

The fitter infers grid size with `int(rows**(1/3)) + 1`, which turns 27 rows into dimension 4 and makes the reshape fail. Recover and validate the integer cube dimension. After changing only that expression in memory, a synthetic 64-row Gaussian fit recovered the input radii and then failed at `create_dataset` because `h5py.File(datafile_name)` opened the file read-only. Make the output destination or write mode explicit. The fitter also hard-codes a different set of `K_T` keys from the HDF5 averager, so those tools need a shared bin convention for end-to-end use. [Grid and file mode](../ebe_scripts/fit_HBT_radii.py#L37) · [Dataset write](../ebe_scripts/fit_HBT_radii.py#L98) · [Averager bin list](../ebe_scripts/average_event_HBT_correlation_function_h5.py#L28)

### 8. Mixed-event analyses lack a policy for missing references and sparse bins

When real mixed events are enabled and their stream yields no events for a nonempty primary batch, HBT and the balance function both use the mixed count as the divisor of a modulo expression. Check the reference batch before sampling and define how an exhausted stream should be handled. Independently, the balance-function output divides by per-bin and total mixed-event counts without a zero-count check. The bundled two-event gzip fixture completed with nonfinite entries in all 21 Δy rows. Keep the raw counts and define how undefined bins are represented and propagated to ensemble outputs. The zero mixed-batch path was identified by source inspection. [HBT mixing](../src/HBT_correlation.cpp#L197) · [Balance-function mixing](../src/BalanceFunction.cpp#L159) · [Δy normalization](../src/BalanceFunction.cpp#L248)

## Verification and reproduction

- Built from the repository root with `cmake -S . -B build/review` and `cmake --build build/review --parallel 2`. Configuring with `-Dunit_tests=ON` stopped because GoogleTest is unavailable. The focused checks used GNU C++ 15.2, NumPy 2.5.4, SciPy 1.18.1, and h5py 3.16.0.
- The Δη₁₃ and decay checks called public C++ functions from small C++11 drivers linked against `libhadronic_afterburner_tools`. For Δη₁₃, load the root `parameters.dat`, set `particle_monval=9999`, `order_max=10`, `n_rap=3`, `rapidity_dis_min=-1`, and `rapidity_dis_max=1`. Use identical `3 × 10` Q-vector grids for all three arguments with `Q_n = 3 exp(i n π/3)` in the middle bin and zeros elsewhere, and call `calculate_three_particle_correlation_deltaeta` with `flag=2` and `flag_ch=1`. For the decay, construct `particle_decay` with the shipped strong table and seed `24680`, look up masses for IDs `113`, `211`, and `-211`, set the mother's `E` to its mass and all three momentum components to zero, then call `perform_two_body_decay`.
- Stage the root `parameters.dat`, `EOS/`, and `unit_tests/test_gzip_reader/particle_samples.gz` in a scratch run directory, with the sample under `results/`. Use the two executable invocations below, substituting an absolute path to the build. For the HDF5 averaging check, load the resulting 27-row, five-column `HBT_correlation_function_KT_*.dat` tables into 20 event groups with multiplicities `1…20`. Provide the six flow datasets named in the averager's `required_files_list` using one row `[0, multiplicity]` for each. Run `python3 ebe_scripts/average_event_HBT_correlation_function_h5.py events.h5 averages` with the script resolved from the checkout and an existing scratch `averages/` directory. To isolate the schema check, resolve only the scalar centrality `min` as `builtins.min` in a scratch copy.
- For the radii script, create a scratch HDF5 group `00-05` with an eight-column dataset named `HBT_correlation_function_KT_0_0.2.dat`, and run `python3 ebe_scripts/fit_HBT_radii.py fit.h5`. A 27-row dataset reaches the reshape failure. For the output-mode check, use a `4 × 4 × 4` mesh with each axis spanning `−0.03` to `0.03 GeV`, populate columns 0–2 with the coordinates, column 6 with the script's Gaussian for `(λ, R_out, R_side, R_long, R_os, R_ol) = (1, 5, 4, 6, 0.1, 0.2)`, and column 7 with `0.02`. Resolve `nq` to the rounded integer cube root in a scratch copy to reach the write failure.

```sh
/absolute/path/to/hadronic_afterburner_tools.e read_in_mode=10 analyze_flow=0 analyze_HBT=1 qnpts=3 randomSeed=12345
/absolute/path/to/hadronic_afterburner_tools.e read_in_mode=10 analyze_flow=0 analyze_balance_function=1 randomSeed=12345
```

The checks used synthetic cases and bundled fixtures. The impact on existing production results was not measured.

## Fixes and follow-up verification

- Reworked the rapidity-differential triplet estimator to integrate the correct harmonic and subtract overlaps according to the particle selections and rapidity bins. Added an ordered-triplet enumeration test for Δη₁₂ and Δη₁₃ with inclusive, same-sign, and opposite-sign selections.
- SMASH header and record reads now report incomplete input and retain independent format versions for primary and mixed streams. A two-event primary v6 stream with a v7 mixed stream produced two events per stream and the expected single same-event pair with a cosine weight of 1.
- HBT now emits a six-column raw companion table while retaining the existing compact and standard table layouts. The averager prefers the companion and also accepts older tables with raw counts. The fitter discovers momentum bins in the ensemble file and masks bins without a usable uncertainty. Verified the new output using the bundled iSS fixture, and ran averaging and fitting on a synthetic Gaussian with a custom momentum bin. The recovered out, side, and longitudinal radii were 5, 4, and 6 fm.
- Percentile scripts use one shared multiplicity-bound implementation. Ties stay together on the peripheral side of the threshold and the last class includes the minimum multiplicity. In the HBT smoke runs, all 20 distinct-multiplicity events were included, and all eight equal-multiplicity events went into the same class.
- Added tests for four-momentum conservation in two- and three-body decays at rest and in motion, and for HBT phase invariance under a rapidity shift.
- HBT and balance-function analyses report an error for an empty external reference batch. Undefined balance-function ratios are explicitly recorded as `NaN`, with raw counts retained in both projections and the two-dimensional output. Verified the C++ output and HDF5 averaging with the bundled fixture.

GoogleTest 1.17 was built from local sources at `/usr/src/googletest`. All 23 tests passed across the three unit executables after refreshing the test parameters and running the flow test through `PerformAnalysis`. Numerical validation used synthetic and bundled samples.
