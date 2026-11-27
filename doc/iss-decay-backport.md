# iSS resonance decay backport

The toolkit now draws decay channels in proportion to their branching ratios over the channels the mother can reach. Four- and five-body decays use Raubold–Lynch phase space, and channels marked with a negative daughter count can open through truncated Breit–Wigner sampling of broad daughters. Exact thresholds and mothers below every reachable channel use the same selection rules as the cascade, so an undecayable mother is retained as a final particle. The implementation uses the toolkit's double-precision masses and existing two- and three-body decay routines. [iSS PR #6](https://github.com/chunshen1987/iSS/pull/6)

The four η₂(1645) → K K̄ π entries in `EOS/pdg-s95pv1.dat` now declare three daughters. The K⁺ K̄⁰ entry uses π⁻ to conserve charge. [iSS PR #7](https://github.com/chunshen1987/iSS/pull/7)

## Verification

- All 25 C++ unit tests passed. The particle-sample tests use the default strong UrQMD table. With seed `12345`, the new tests check four-body and below-pole channels, four-momentum and quantum numbers, and the exact-threshold and terminal-resonance cases.
- A C++11 driver linked against `libhadronic_afterburner_tools` sampled 102,200 cascades across the bundled strong and weak UrQMD, s95pv1, and SMASH tables. Each run pointed a scratch `EOS/pdg.dat` at the selected table. The driver used seed `97655`, mother masses at their poles, and longitudinal momentum of `0.37 GeV`. It ran 100 cascades per unstable table entry and 10,000 for η₂(1645) where present. It verified cascade termination, on-shell daughters, four-momentum conservation to `1e-10 GeV`, and conservation of charge and baryon number. All four corrected η₂ channels were observed.
- A scratch table supplied a 1 GeV mother at rest with a five-photon channel of weight `0.2` and a closed channel to two 0.6 GeV daughters of weight `0.3`. With the same seed, all 10,000 draws produced five photons. The mean `2E/M` of the first photon was `0.399356`, compared with the massless five-body phase-space expectation of `0.4`.

## November 26 optimization branch

The two-, three-, and n-body routines now calculate the rest-frame mass from the mother's stored four-momentum and use a boost expressed directly in terms of that four-momentum. The same rest-frame mass sets the time dilation for the decay point. A channel selected with the listed mass retains that mass when rounding of the stored four-vector would close the channel. The high-boost regression uses float-rounded mother momenta at Lorentz factors 3, 30, and 300 and checks four-momentum conservation and daughter mass shells. [iSS fast-resonance fix](https://github.com/chunshen1987/iSS/commit/e371aa06d473ecab78db0e34cc001c9b774d8e44)

The toolkit reuses decay scratch space and can process independent events with OpenMP. Set `resonance_decay_parallel=1` and use `OMP_NUM_THREADS` to choose the number of workers. This mode derives one random stream per event from the batch seed, so a fixed seed and event buffer size give identical results across worker counts. The default uses the existing sequential stream. Builds without OpenMP can run the event-stream mode sequentially. [iSS parallel-decay change](https://github.com/chunshen1987/iSS/commit/74934891123c84726fd485ac29ded1eff3b9f199)

All 27 C++ tests passed in both the OpenMP build and a build configured with `-DCMAKE_DISABLE_FIND_PACKAGE_OpenMP=TRUE`. With seed `4`, the parallel regression decayed eight events of 32 f₂(1270) mothers and compared every daughter's identity, four-momentum, and decay time at one and four workers. The results were identical.
