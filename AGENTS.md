# AGENTS.md

## Scope and code map

This toolkit analyzes final-state hadron samples from high-energy nuclear collisions. Its observables include identified and charged-particle spectra, anisotropic flow and multiparticle correlations, Hanbury Brown–Twiss (HBT) femtoscopy, charge balance functions, and event-by-event yields.

- `src/main.cpp` and `src/Analysis.cpp` load `parameters.dat`, apply command-line `key=value` overrides, and dispatch analysis modules. The executable uses `results/` for event input and analysis output, relative to the working directory.
- `src/particleSamples.cpp` handles input formats, event buffering, particle selection, and kinematics. It works with `src/particle_decay.cpp` and `EOS/` for resonance decays and feed-down. Record the actual targets of `EOS/pdg.dat` and `EOS/pdg_weak.dat` when comparing decay results.
- `src/single_particleSpectra.cpp`, `src/HBT_correlation.cpp`, `src/BalanceFunction.cpp`, and `src/particle_yield_distribution.cpp` implement the observable calculations. `ebe_scripts/` contains format converters, event collection scripts, and ensemble post-processing, including HDF5 aggregation and HBT fits.

## Physics and numerical conventions

- For changes to an observable, state the estimator, acceptance cuts, reference particles or subevents, and normalization in the change summary. Cite the relevant equation or reference when implementing a published definition. Keep the existing defaults unless the task calls for a physics change, and quantify its effect on representative output when feasible.
- Track units at reader and observable boundaries. Momentum, energy, and mass are in GeV. Spacetime coordinates are in fm with `c = 1`. Check conversions for each input format and use the existing `hbarC` or `hbarc` convention in the affected calculation. The HBT phase `q · Δx / (ℏc)` must be dimensionless.
- Distinguish momentum rapidity `y`, pseudorapidity `η`, and spacetime rapidity. Follow the observable's `rap_type`, rapidity shift, azimuth convention, and bin-edge inclusivity. For HBT, verify the signs of `q = p1 - p2`, the definition `K = (p1 + p2)/2`, and the frame and out–side–long axes before altering pair kinematics.
- Preserve particle ID and charge-selection semantics, including generator-to-PDG mappings, antiparticles, aggregate particle selectors, and resonance feed-down flags. For decay changes, check kinematic thresholds and four-momentum and charge conservation where applicable.
- Preserve physical event boundaries and distinguish oversamples of the same collision from independent collisions. For correlations, check self-pair subtraction, ordered versus unordered combinations, multiplicity weights, reference acceptance, and same-event versus mixed-event normalization. In ensemble scripts, preserve the centrality estimator and the independent event used for uncertainty resampling. A mean of per-event ratios can differ from a ratio of summed numerators and denominators.
- Use double precision for kinematics and sensitive accumulations. Check empty bins, low multiplicities, zero denominators, and logarithm and square-root domains in the code being changed. Preserve the established treatment of undefined bins, or explain and validate any change. Compare floating-point results with tolerances tied to numerical precision or Monte Carlo uncertainty.

## Implementation and compatibility

- Keep C++ code compatible with C++11 and follow `.clang-format` for edited code. Reuse the existing readers, `ParameterReader`, and `RandomUtil::Random` where appropriate. Keep event buffering bounded, and check event or pair counts and runtime if changing hot loops or buffer handling.
- Treat `read_in_mode`, input filenames and layouts, parameter names, text output columns, and HDF5 dataset names as interfaces. If a change touches one, trace its producers and consumers in `src/` and `ebe_scripts/`. Check the relevant text, gzip, or binary format using a matching fixture. For binary changes, verify record layout, field widths, and behavior at EOF and on short reads.
- Keep generated event files, build products, and large output datasets out of commits. Run executable and post-processing checks in a scratch run directory with copied inputs because output may be written alongside inputs.

## Build and validate

CMake requires a C++11 compiler and zlib. From the repository root, use an out-of-source build. `compile.sh` also reformats sources, installs into the source tree, and removes its build directory.

```sh
cmake -S . -B build
cmake --build build --parallel 2
```

The executable is `build/src/hadronic_afterburner_tools.e`. For a smoke run, stage `parameters.dat`, `EOS/`, and a `results/` directory containing the input files for the chosen `read_in_mode`. Invoke the executable by absolute path from that run directory.

GoogleTest is required for `-Dunit_tests=ON`. The three test programs use fixture paths relative to the working directory. Run them directly from a copy of `unit_tests/`, with `EOS/` alongside it so the relative link resolves. This keeps test output out of the source fixtures.

```sh
cmake -S . -B build/tests -Dunit_tests=ON
cmake --build build/tests --parallel 2
repo_root=$PWD
test_run=$(mktemp -d)
cp -R EOS unit_tests "$test_run/"
(
    cd "$test_run/unit_tests" || exit 1
    "$repo_root/build/tests/unit_tests/particleSamples_unittest.e" &&
    "$repo_root/build/tests/unit_tests/HBT_unittest.e" &&
    "$repo_root/build/tests/unit_tests/Analysis_unittest.e"
)
```

- For physics changes, use a small sample with a known result where feasible. Examples include hand-counted particles at acceptance edges, explicit particle-tuple enumeration for a Q-vector correlator, and a known pair phase or longitudinal boost for HBT. Check the relevant output values and normalization in addition to successful execution.
- Use a fixed `randomSeed` for stochastic regression checks. Record the input format and sample, parameters and overrides, and seed so the comparison can be repeated. For a post-processing change, exercise the affected script on a small representative dataset and check the output schema and numerical result. Report any validation blocked by missing dependencies or fixtures.

## Investigation reports

For each investigation task, save a Markdown summary report in `doc/` and link to it in the final response. Use a descriptive filename and record the question, key findings, supporting evidence, checks performed, and any unresolved limitations. Include the inputs, parameters, and commands needed to reproduce numerical findings when applicable.
