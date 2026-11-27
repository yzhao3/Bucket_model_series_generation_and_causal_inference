# Hydrometeorological Bucket Model: Time-Series Generation and L–K Causality Analysis

This repository provides MATLAB code for generating synthetic time series from a hydrometeorological bucket model and estimating multivariate Liang–Kleeman (L–K) information flow to reconstruct the embedded causal relationships.

## 1. Files

| File | Description |
| --- | --- |
| `example.m` | Example driver: sets simulation parameters, generates time series, adds optional observational noise, and estimates directional information flows within a four-variable multivariate system. |
| `bucket_model_construction.m` | Generates rainfall (`R`), soil water content (`S`), lateral interflow (`I`), and runoff (`Q`) time series. |
| `multi_causality_est.m` | Estimates information flow from the **second** input series to the **first** within the multi-dimensional matrix, with uncertainty thresholds at 90%, 95%, and 99% confidence levels. |

## 2. Software requirements

- MATLAB. The associated manuscript reports use of **MATLAB R2021b**.
- **Econometrics Toolbox** for `dtmc` and `simulate` used to generate rainfall occurrence.
- **Statistics and Machine Learning Toolbox** for `betarnd` used to generate rainfall amounts.
- The three `.m` files should be in the same working directory, or their location should be on the MATLAB path.

No external observational or reanalysis dataset is required for **this synthetic bucket-model example**. The driver generates its own time series. The `example.m` script contains no file-reading or file-writing commands.

## 3. Model inputs and generated time series

The bucket model generates four column-vector time series, each containing `num_steps` values:

- `R`: rainfall forcing;
- `S`: soil water content;
- `I`: lateral interflow;
- `Q`: runoff.

Model settings are defined in `bucket_model_construction.m`, including maximum storage `S_max = 80`, initial storage `S_0 = 40`, depletion coefficient `Ks = 2.3`, interflow threshold `delta = 10`, exponent `xi = 0.6`, and process-noise lag-one correlation parameter `r = 0.8`.

Rainfall occurrence is sampled from a two-state Markov chain with transition matrix `[0.8 0.2; 0.2 0.8]`; rainfall amounts use a scaled beta distribution (`alpha = 9`, `beta = 6`).

The driver first creates a simulation with the process-noise powers set to zero to estimate reference signal powers. It then runs simulations with process noise scaled by `pn_SNR`. Optional independent additive Gaussian observational noise is controlled by `on_switch` and `on_SNR`.

## 4. Quick start

1. Place `example.m`, `bucket_model_construction.m`, and `multi_causality_est.m` in one directory.
2. Start MATLAB, navigate to that directory, and run:

   ```matlab
   run('example.m')
   ```

3. After completion, inspect the variables in the MATLAB workspace:

   ```matlab
   labels = {'R','S','I','Q'};
   disp(multi_T21)       % Estimated signed information-flow matrix
   disp(multi_err)       % 99% uncertainty thresholds
   disp(sig_T > 0)       % Binary significance matrix (99% level)
   ```

**Direction convention:** In the supplied `example.m`, `multi_T21(i,j)` is the estimated flow **from variable `i` to variable `j`**. Rows are sources (transmitters) and columns are targets (receivers), in the order `[R, S, I, Q]`. For instance, `multi_T21(1,2)` means `R -> S`. This matches the transmitter-row / receiver-column layout of the manuscript's information-flow table.

The matrix diagonal remains zero because self-links are not estimated.

## 5. Parameters in `example.m`

| Parameter | Supplied default | Meaning |
| --- | --- | --- |
| `num_steps` | `10000` | Number of simulated time steps in one realization. |
| `nsim` | `100` | Number of simulated realizations. |
| `dt` | `1` | Time interval specified in the driver; the estimator also uses a hard-coded internal `dt = 1`. |
| `IF_length` | `3000` | Number of final time steps retained for causality estimation. |
| `pn_SNR` | `10000` | Process-noise signal-to-noise ratio. |
| `on_SNR` | `10` | Observational-noise signal-to-noise ratio when noise is enabled. |
| `on_switch` | `1` | `1`: add observational noise; `0`: do not add observational noise. |

The example estimates information flow on the last `IF_length` samples of each generated realization, using a forward-difference time advance of `np = 1` in the estimator.

### Example configuration for the manuscript's low-process-noise benchmark

In `example.m`, set:

```matlab
IF_length = 3000;
pn_SNR = 10000;
on_switch = 0;
nsim = 100;
```

**Important:** Matching these parameter values does **not** automatically reproduce the manuscript's mean information-flow table. The supplied driver currently overwrites `multi_T21`, `multi_err`, and `sig_T` on each simulation and retains only the **final** realization. To reproduce across-realization means and uncertainty summaries, the results need to be saved separately for each realization and summarized after the loop.

## 6. Computation workflow

1. **Estimate reference signal powers.** Call `bucket_model_construction` with zero process-noise powers; calculate `mean(R.^2)`, `mean(S.^2)`, `mean(I.^2)`, and `mean(Q.^2)`.
2. **Generate a noisy realization.** Call the bucket-model function using the selected `pn_SNR` and reference powers.
3. **Optionally add observational noise.** When `on_switch = 1`, add independent Gaussian errors scaled by the corresponding variable's signal power and `on_SNR`.
4. **Extract the analysis window.** Form an `IF_length × 4` matrix whose columns are `[R, S, I, Q]`, using the final `IF_length` time steps.
5. **Estimate directional flows.** For each source `i` and target `j` with `i ~= j`, rearrange the columns as `[target(j), source(i), remaining variables]`, as done in the driver, and call `multi_causality_est`. The estimator computes flow from its second input column to its first, so the result stored in `multi_T21(i,j)` represents `i -> j`.
6. **Evaluate significance.** Store the 99% uncertainty threshold in `multi_err(i,j)`. The driver calculates `sig_T = abs(multi_T21) - abs(multi_err)`; positive entries are treated as significant at the selected 99% level.
7. **Compare networks.** Use significant directed links to summarize the inferred network and compare it with the bucket model's expected causal relationships.

## 7. Outputs and interpretation

The supplied `example.m` leaves the following in the MATLAB workspace for its **last** iteration:

| Variable | Dimensions | Description |
| --- | --- | --- |
| `R`, `S`, `I`, `Q` | `num_steps × 1` | Synthetic hydrometeorological time series. |
| `varout` | `num_steps × 4` | Time series with optional observational noise. |
| `varif` | `IF_length × 4` | Time-series segment used for information-flow analysis. |
| `multi_T21` | `4 × 4` | Signed directed information-flow estimates. |
| `multi_err` | `4 × 4` | Corresponding 99% uncertainty thresholds. |
| `sig_T` | `4 × 4` | `abs(multi_T21) - abs(multi_err)`, not a p-value. |

The output is **not automatically plotted or saved to disk**. Numerical values vary between runs because the time series are randomly generated unless a seed is specified.

For the bucket-model configuration described in the associated manuscript, expected directed links are `R -> S`, `R -> Q`, `S -> I`, `S -> Q`, `I -> S`, `I -> Q`, and `Q -> I`. Recovery of these links is an empirical validation target; it should not be treated as guaranteed for every short or noisy realization.

## 9. References

- Ombadi, M., et al. (2020), *Evaluation of Methods for Causal Discovery in Hydrometeorological Systems*, Water Resources Research, 56, e2020WR027251. https://doi.org/10.1029/2020WR027251
- Liang, X. S. (2021), *Normalized Multivariate Time Series Causality Analysis and Causal Graph Reconstruction*, Entropy, 23, 679. https://doi.org/10.3390/e23060679
- Liang, X. S. (2016), *Information flow and causality as rigorous notions ab initio*, Physical Review E, 94, 052201. https://doi.org/10.1103/PhysRevE.94.052201
