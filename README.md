hadronic_afterburner_toolkit
===========================

## This program deals with hadronic afterburner output files and perform various analysis on hadronic observables

### Zlib is required to read in zipped UrQMD files

### Charged and Identified particle spectra and flow coefficients

Single particle spectra for charged hadrons and identified particles are 
collected.
The anisotropic flow coefficients from two particle correlation are also 
analyzed. As well as event-plane correlation, vn_distribution, rn ratio
for flow factorizatin breaking.

### Particle pair HBT correlation from Monte-Carlo samples of emitted particles

The HBT correlation between identified particle pair is originated from the quantum statistic of the final wave function. The Monte-Carlo samples of the emitted particles do not have this correlation built in. So we need to add it by hand. 

The HBT correlation function can be approximatly computed as,

C(q, K) - 1 = \frac{\int d^4 x \int d^4 y s(x, p1) s(y, p2) cos(q \cdot (x - y)}{\int d^4 x \int d^4 y s(x, p1) s(y, p2)},

where K = 0.5(p1 + p2) and q = p1 - p2. 

In the program, we introduce oversampling parameter to increase statistics by first grouping multiply events together to a single big event and then form the pairs. The number of the pairs will increase by oversampling factor.

Pairs from mixed event are form by first rotating the mixture event by a random angle which eliminate the unwanted correlations.

The output correlation function is in 3-d q space (q_out, q_side, q_long).

For HDF5 ensemble averaging, leave `hbt_raw_output=1` and collect the
`HBT_raw_correlation_function_KT_*.dat` files alongside the existing HBT
output and flow files used for centrality selection. The raw tables store
three nominal q-bin coordinates in GeV, same-event pair count, sum of the
cosine weights, and mixed-event pair count. The raw files are available in
both standard and compact output modes.

Run `ebe_scripts/average_event_HBT_correlation_function_h5.py events.h5 output_dir`
on the combined event file. Then run `ebe_scripts/fit_HBT_radii.py` on
the resulting ensemble HDF5 file. The fitter writes radii into that file
using its existing momentum-bin names. Regenerate old two- or five-column
files for ensemble normalization if the raw pair counts were not saved.

### Event-by-event particle yield distribution

Particle yield distrbution from a given pT and rapidity cut can be collected from UrQMD outputs

### Parallel resonance feed-down

With resonance feed-down enabled, set `resonance_decay_parallel=1` to
assign independent random streams to events. CMake enables OpenMP when
available, and `OMP_NUM_THREADS` controls the worker count. For a fixed
`randomSeed` and event buffer size, the results in this mode are independent
of the worker count. The default `resonance_decay_parallel=0` uses a single
sequential random stream.

