# Parameters and CSV output guide

This document describes the public experiment scripts only. The original research report uses `ω` for the dynamic rewiring threshold; the code uses `TH`.

## Parameters

| Code name | Phase(s) | Definition |
|---|---|---|
| `N` | I.S, I.D, II.S | Square-lattice side length, giving `N*N` agents |
| `T` | I.S, I.D, II.S | Simulation duration in discrete time steps |
| `p` | I.S, I.D, II.S | Initial network rewiring probability |
| `A` | I.S, I.D, II.S | Fraction of information redistributed in threshold-triggered activation |
| `f` | I.D | Rewiring frequency setting |
| `TH` | I.D | Edge-use threshold for rewiring (report: `ω`) |
| `B` | II.S | Secondary probabilistic contagion sensitivity |
| `s_min` | All | Lower cutoff used for avalanche tail fitting |

## Common summary columns

| CSV column | Explanation |
|---|---|
| `n_avalanches` | Number of avalanches retained by the simulation summary |
| `n_tail` ($\eta_{tail}$) | Number of avalanches above the chosen lower fitting cutoff |
| `tau_hat` ($\hat{\tau}$) | Estimated avalanche-size distribution tail exponent |
| `ks` ($D_{KS}$) | Kolmogorov–Smirnov distance between empirical and fitted distributions |
| `size_p99` ($S_{0.99}$) | 99th percentile avalanche size |
| `size_max` ($S_{max}$) | Maximum observed avalanche size |
| `frac_system_spanning` ($\phi_{span}$) | Fraction of avalanches activating at least 90% of the agents |

Phase I.D additionally reports clustering, transitivity, assortativity, modularity, rewiring, and edge-persistence diagnostics. Phase II.S additionally reports shares, ratios, prevalence, and correlations of secondary (`P`) activation. Consult the exact CSV header and source functions for each specific derived quantity; this glossary does not replace a code-level definition.


## Symbols in the selected later-phase figures

The Phase II.D, III.S and III.D images are included to illustrate the wider study; these labels refer to the mathematical notation in the research report, not necessarily columns present in the published datasets.

| Figure notation | Explanation |
|---|---|
| $\hat{\tau}$ | Fitted avalanche-size power-law tail exponent |
| $D_{KS}$ | Kolmogorov–Smirnov distance between the empirical and fitted distributions; smaller indicates closer agreement |
| $\phi_{\mathrm{span}}$ | Fraction of avalanches reaching at least 90% of agents |
| $S_{0.99}$ | 99th percentile of avalanche sizes |
| $S_{\max}$ | Largest recorded avalanche size |
| $\langle c\rangle_T$, $\langle c\rangle_F$ | Mean local clustering coefficient computed separately within the traditional-agent (T) and FOMO-agent (F) layers |
| $f$ | Rewiring frequency |
| $\omega$ | Rewiring threshold |
| $p_s$ | Probability of cross-layer connections |
| $A_s$ | Strength of cross-layer redistribution |
| $\pi_F$ | Fraction of agents assigned to the F-layer |

Each heatmap's colour bar measures the quantity named above that panel. Colours should be interpreted against that panel's own scale.
