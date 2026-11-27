# Scientific code review of revision 9ce7e81

Reviewed the recent decay backport and the reader, selection, and observable paths that consume its output. Findings below concern behavior still present at this revision. All 27 existing C++ tests passed in the OpenMP build.

## High-priority findings

### 1. Net-particle yields use the wrong rapidity for antiparticles

The antiparticle loop reads rapidity and pseudorapidity from the particle list, then reads transverse momentum from the antiparticle list. A synthetic event with a positive pion at $y=1$ and a negative pion at $y=0$, using the default $-0.5<y<0.5$ cut, produced $N_+-N_-=0$ instead of $-1$. If the antiparticle list is longer, the same code can index beyond the particle list. Use the same antiparticle record for every cut. [Antiparticle loop](../src/particle_yield_distribution.cpp#L88) · [Particle accessors](../src/particleSamples.h#L244)

### 2. HBT treats the read buffer as one collision for same-event pairs

The numerator pools all event records in a buffer into one particle list before pairing. For a file with two independent events and one pion per event, the raw same-event count was 1 with both events in the same buffer and 0 with `event_buffer_size=1`. Each independent event has zero same-event pairs. Pooling oversamples can be meaningful when they share a physical collision, but the buffer does not encode that relationship. Preserve event boundaries for independent samples and make oversample grouping explicit. [Event grouping](../src/HBT_correlation.cpp#L191) · [Pooled pair loop](../src/HBT_correlation.cpp#L257)

### 3. HBT applies the documented pair-rapidity cut to each track

`HBTrap_min` and `HBTrap_max` are described as pair-rapidity bounds, but both paths discard tracks before forming pairs. A pair with equal and opposite longitudinal momenta has pair rapidity zero and should pass the default bounds regardless of whether one or both tracks lie outside them individually. In a synthetic event, pions at approximately $y=+0.60$ and $y=-0.60$ produced zero same-event pairs, while pions at approximately $y=+0.45$ and $y=-0.45$ produced one. Test the pair kinematics in both numerator and denominator, and apply single-track acceptance independently if needed. [Parameter definitions](../parameters.dat#L114) · [Same-event filter](../src/HBT_correlation.cpp#L261) · [Mixed-event filter](../src/HBT_correlation.cpp#L474)

### 4. Balance functions discard distinct particles at equal rapidity

Both pair loops skip every pair with `abs(delta_y) < 1e-10`, even when the tracks are distinct and have different charges. A synthetic unlike-sign pion pair at $y=0$ gave zero opposite-sign raw counts. Shifting one of those tracks to $y=0.05$ gave two entries, one for each charge orientation. This creates a central-bin deficit for sources or inputs with repeated rapidities. Exclude self-pairs by track identity in overlapping same-event selections. Preserve distinct equal-rapidity pairs in both paths. [Same-event veto](../src/BalanceFunction.cpp#L158) · [Mixed-event veto](../src/BalanceFunction.cpp#L202)

## Additional findings

### 5. The phi multiplicity histogram does not sum to one

The histogram divides every `P(N)` bin by the $K^+K^-$ branching fraction of 0.489, even in iSS binary mode, which disables phi reconstruction. A one-event input containing one accepted, explicitly stored phi wrote `P(1)=2.04498978`. The histogram of observed multiplicities should sum to one. Recovering a parent multiplicity distribution from reconstructed daughters requires an efficiency model, since a branching fraction changes the probability of migrating between multiplicity bins. [Reconstruction flag](../src/particleSamples.cpp#L440) · [Normalization](../src/particle_yield_distribution.cpp#L154)

### 6. Two feed-down controls do not change the active calculation

`select_resonances_flag` loads and prints a list, but the cascade never consults it. Sixty single-ρ⁰ events gave byte-identical charged-pion histograms with the selector enabled and disabled, even though the configured list contains only PDG 213. Separately, the Σ⁰ → Λγ switch populates the parent list, but its decay function has no caller. This second path was identified by source inspection. Connect the settings to the active cascade or filtering path and test selected and unselected parents. [Selector configuration](../parameters.dat#L31) · [List load](../src/particleSamples.cpp#L376) · [Cascade](../src/particleSamples.cpp#L2111) · [Sigma parent selection](../src/particleSamples.cpp#L1311) · [Sigma decay function](../src/particleSamples.cpp#L2183)

### 7. Decay threshold handling can produce off-shell or nonfinite daughters

When the stored mother four-vector is below threshold but its listed mass opens the channel, the fallback uses the listed mass to create rest-frame daughters and then uses the inconsistent stored four-vector for the boost. This conserves the stored total by allowing individual daughters to leave their mass shells. With a test mother at rest with listed mass 0.30 GeV, stored energy 0.25 GeV, and two 0.138 GeV daughters, the calculation gave $E^2-\mathbf p^2=0.012169\;\mathrm{GeV}^2$ instead of $0.019044\;\mathrm{GeV}^2$. A separate direct call to the public three-body routine at exact threshold (0.375 GeV into three 0.125 GeV daughters) gave two nonfinite daughters after the angular calculation divided by zero. Normal table selection excludes exact narrow-daughter thresholds, but the public routine accepts equality. Make the policy for incompatible input four-vectors explicit and handle zero available phase space before angular sampling. [Rest-mass fallback and boost](../src/particle_decay.cpp#L404) · [Three-body threshold and angular sampling](../src/particle_decay.cpp#L616)

### 8. HBT fitting constrains both cross terms to nonnegative values

The signs of `R_os` and `R_ol` disappear because the parameters are squared. A Gaussian with $R_{ol}^2=-4\;\mathrm{fm}^2$ gives correlations 0.29376 and 0.20295 at $(q_o,q_s,q_l)=(0.03,0,+0.03)$ and $(0.03,0,-0.03)$ GeV, respectively, for $R_o=5\;\mathrm{fm}$, $R_l=6\;\mathrm{fm}$, and unit intercept. The current parametrization can only produce the reverse ordering or equality. This can bias the diagonal radii in bins where a signed cross term is present. Use signed cross coefficients with a positive-definite width matrix and state their units in the output. [Fit model](../ebe_scripts/fit_HBT_radii.py#L13)

### 9. Flow analyses leak their spectra objects and retain the final sample buffers

The flow methods allocate spectra with `new` and clear only the vector of pointers at the end. Each spectra object also keeps a shared pointer to the particle sample. Repeating analyses in a long-lived process retains both the spectra and their final event buffers. Identified by ownership inspection. Use owning smart pointers for the spectra objects. [Allocation](../src/Analysis.cpp#L123) · [Cleanup](../src/Analysis.cpp#L321) · [Retained sample](../src/single_particleSpectra.cpp#L259)

## Reproduction notes

- The CLI checks used `build/fixes-tests/src/hadronic_afterburner_tools.e` in scratch directories with the root `parameters.dat` and `EOS/`. Each `results/particle_list.bin` contained a native 32-bit particle count per event and `struct.pack("=i9f", pdg, m, t, x, y, z, E, px, py, pz)` per particle. All positions and times were zero. Unless stated otherwise, the pions had $m=0.138$ GeV, positive $p_x=0.23$ or $0.24$ GeV, and $p_y=0$. Set $E=m_T\cosh y$ and $p_z=m_T\sinh y$. Run with `read_in_mode=9 analyze_flow=0 randomSeed=12345`.
- For HBT, use `analyze_HBT=1 qnpts=11 q_min=-1 q_max=1`, then sum column 3 (zero-indexed) in the raw tables. For the buffer check, the two pion momenta were $p_x=0.23,0.24$ GeV, one per event. For the rapidity check the same momenta were placed in one event with opposite rapidities. Compare `event_buffer_size=1` with the default.
- For the balance function, use `analyze_balance_function=1` and sum column 3 of `Balance_function_9998_-9998_Delta_y.dat`. For the net yield, use `analyze_ebe_yield=1 net_particle_flag=1` and inspect the occupied bin in `particle_0211_yield_distribution_y.dat`. For the phi yield, use `analyze_ebe_yield=1 particle_monval=333` with $m=1.019$ GeV. For the resonance-selector comparison, use PDG 113 at $m=0.769$ GeV in each of 60 events and run `analyze_ebe_yield=1 resonance_feed_down_flag=1 rap_min=-3 rap_max=3` with both values of `select_resonances_flag`.
- The decay checks used a C++11 driver linked against `libhadronic_afterburner_tools` with the strong UrQMD table and seed `12345`. Call `perform_two_body_decay` with mother PDG 113 and `perform_three_body_decay` with mother PDG 223 using the masses and energies above. The odd cross-term check evaluated `gaussian_3d` from `fit_HBT_radii.py` directly with the two reflected momentum points.

Validation used synthetic events, the existing unit tests, and source inspection. No production results were recomputed.

## Fixes and verification

- Net yields now apply every acceptance cut to the same track. The original rapidity case and a pseudorapidity case containing only an antiparticle both gave `P(-1)=1`. The stored-φ test gave `P(1)=1`.
- HBT now accumulates pairs separately for each collision, with explicit grouping of consecutive oversamples. A four-event fixture with one pion per record gave same-event counts of 0, 2, and 3 for group widths 1, 2, and 3. The same and mixed counts were identical across buffer sizes 1, 2, and 100. Without a separate mixed stream, each collision is mixed against the preceding collision after a random azimuthal rotation. A separate-stream test with tracks outside the single-track rapidity bounds accepted one same-event and one mixed-event pair inside the pair-rapidity bounds.
- The balance-function test used two identical positive pions and one negative pion at equal rapidity. It counted four unlike-sign and two like-sign ordered pairs. Only a track paired with itself was removed in the same-event calculation.
- With selective feed-down, a PDG 213 parent in the configured list produced a positive pion and an unselected PDG 113 parent was kept. The dedicated Sigma switch produced a Λ or anti-Λ from the corresponding signed parent when enabled.
- Exact-threshold three-body daughters have zero rest-frame momentum. A rounding-sized discrepancy at threshold is repaired by recomputing the parent energy from its listed mass and stored momentum. The large inconsistency from the review raises a diagnostic instead of emitting off-shell daughters.
- A synthetic $13^3$ Gaussian grid with side, out, and long radii of 4, 5, and 6 fm recovered signed out–side and out–long coefficients of +2 and −4 fm² in all five fit intervals. The HDF5 header now labels these columns in fm².
- All 29 C++ tests passed. The focused `Analysis.FlowAnalysis` test also passed under AddressSanitizer with leak detection enabled after the spectra objects gained owning pointers.
