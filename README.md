# Robust Spectral Parametrization 
% Author: Ying Wang, Min Li
% Create Time: 2025
% Copyright (c): 2020-2025 Ying Wang, yingwangrigel@gmail.com,
%                Min Li, minli.231314@gmail.com
% Joint China-Cuba LAB, UESTC, Hangzhou Dianzi Univerisity 
% License: GNU General Public License v3.0 (see LICENSE file)

RobustSpecPara is a MATLAB toolbox for robustly decomposing EEG power spectra into aperiodic (``\xi``) and oscillatory (alpha-band) constituents. It initialize non-linear fits, and export parameter tables for downstream parameterization normative EEG (https://github.com/LMNonlinear/AP-qEEG).

## Highlights
- Implements Lorentzian aperiodic and Gaussian oscillatory kernels that are jointly optimized with Levenberg-Marquardt for stable estimates of the ``\xi`` slope, offsets, peak frequency, bandwidth, and power.
- Provides scripts that process HarMNqEEG cross-spectra and export channel-wise tables of alpha and ``\xi`` parameters for normative analyses.
- Ships quality metrics (e.g., weighted R^2) and helpers for smoothing, bounds handling, and table manipulation.


## Methodological Context

![RobustSpecPara](./image/readme/RobustSpecPara.png)

### Robust EEG Spectral Parametrization
RobustSpecPara implements a principled approach to decompose EEG power spectra into two fundamental components:

1. **Aperiodic (ξ) Component** - The 1/f background activity
   - Modelled using a **Lorentzian kernel** with parameters:
     - **χ (chi)**: Exponent controlling the slope of the aperiodic background
     - **h (height)**: Amplitude of the aperiodic component
     - **μ (mu)**: Knee frequency (transition point)
     - **Baseline offsets**: Intercept terms for flexible baseline fitting
   - Captures the broadband spectral characteristics independent of oscillatory activity
   - Robust to variations in recording conditions and electrode impedance

2. **Oscillatory (α) Component** - Periodic activity (e.g., alpha rhythm)
   - Modelled using **Gaussian kernels** with parameters:
     - **μ (mu)**: Peak frequency (centre of the oscillation)
     - **σ (sigma)**: Bandwidth (width of the peak)
     - **h (height)**: Amplitude of the oscillatory component
   - Captures narrow-band oscillations superimposed on the aperiodic background
   - Multiple peaks can be detected across frequency bands (delta, theta, alpha, beta)

3. **Joint Optimization**
   - Both components are fitted **simultaneously** to the linear-scale PSD
   - Uses **Levenberg-Marquardt** non-linear optimization for stable parameter estimates
   - Avoids multiplicative assumptions made by log-domain approaches (e.g., FOOOF)
   - Elastic constraint regularization prevents overfitting and ensures physiologically plausible solutions


## Requirements
- MATLAB R2021a or newer (tested with >= R2022b).
- Signal Processing Toolbox (for PSD utilities) and Parallel Computing Toolbox (scripts use `parfor`).
- Access to HarMNqEEG spectral tables or your own PSD matrices.

## Installation & Setup
1. Clone or download `RobustXiAlpha` into your MATLAB workspace.
2. Add the toolbox folders to the MATLAB path:
   ```matlab
   addpath(genpath(fullfile(pwd, 'utility')));
   addpath(genpath(fullfile(pwd, 'external')));
   addpath(genpath(fullfile(pwd, 'example')));
   ```
3. If you plan to use parallel execution, confirm that a MATLAB parallel pool can be opened (`parpool`, or the provided `startmatlabpool` wrapper).

## Quick Start

### Option 1: Simple Spectral Decomposition (Recommended for New Users)
Execute `main_robust_s_para.m` for a straightforward example:
```matlab
run('example/main_robust_s_para.m')
```
This script:
- Loads a sample power spectrum from `data/Spec_example.csv` (47 frequency bins, log-scale).
- Converts log-scale PSD to linear scale and applies Nadaraya-Watson smoothing.
- Fits a single BiXiAlpha model decomposing the spectrum into:
  - **Aperiodic (ξ) component**: Lorentzian kernel capturing the 1/f background.
  - **Oscillatory (α) component**: Gaussian kernel capturing the alpha-band peak.
- Saves fitted parameters to `results/robust_s_para/robust_s_para_analysis/xialpha_result.mat`.
- Supports debug mode (single channel) or parallel processing (all channels).

## Input Specification
The main entry point expects:
- **`f`** - Frequency vector (Hz) sampled uniformly, typically 1–47 Hz for HarMNqEEG or custom ranges.
- **`s`** - Power spectral density values (linear scale, not log). Columns represent channels; rows represent frequency bins.
  - Input data should be in **linear scale** (not dB or log). If your data is log-scale, convert using `s_linear = exp(s_log)`.
- **`bs`** *(optional)* - Bispectral terms for higher-order spectral analysis; pass `[]` if not available.
- **`ks`, `kh`** - Taper parameters for spectral smoothing:
  - `ks` (default 4): Number of harmonics for the aperiodic (ξ) component.
  - `kh` (default 8): Number of harmonics for the oscillatory (α) component.
  - Higher values increase model complexity; typical range is 3–10.
- **`Fs`** - Sampling rate (Hz) used for normalisation; typically `2 × max(f)`.

### Data Preparation
- **Log-to-linear conversion**: If your input is log-scale PSD (common in qEEG), convert using `s = exp(s_log)`.
- **Smoothing**: Apply Nadaraya-Watson or other kernel smoothing to reduce noise (see `psd_smooth_nw` helper).
- **Metadata**: Age, sex, country, device, and dataset information are attached to each `BiXiAlpha` instance via the `.info` struct for downstream normative modelling.

## Outputs
For every channel, `BiXiAlpha` stores comprehensive fitted parameters and quality metrics:

### Aperiodic (ξ) Component Parameters
- **`para_xi.kernel.para`** - Lorentzian parameters:
  - **χ (chi)**: Exponent of the aperiodic slope (typically 0.5–2.5)
  - **h (height)**: Amplitude of the aperiodic component
  - **μ (mu)**: Knee frequency (transition point)
  - **Baseline offsets**: Intercept terms for flexible baseline fitting

### Oscillatory (α) Component Parameters
- **`para_alpha.kernel.para`** - Gaussian peak parameters (one or more peaks):
  - **μ (mu)**: Peak frequency (Hz)
  - **σ (sigma)**: Bandwidth (Hz)
  - **h (height)**: Amplitude of the oscillatory component

### Quality-of-Fit Metrics
- **`r2`** - Coefficient of determination (plain R²)
- **`r2w`** - Variance-weighted R² (recommended when observation variances are known)
- **`r2_s`** - Spectral R² (alternative metric)
- **`loss`** - Final optimization loss value
- **`regularization`** - Type of regularization applied (e.g., `elastic_constraint_mu_sigma_chi_xi`)
- **`lambda`** - Regularization weights for different parameters

### Saved Results
Saved `.mat` files preserve arrays of `BiXiAlpha` objects that can be reloaded for:


### Key Strengths
1. **Simultaneous Decomposition** - Aperiodic and oscillatory components are fitted jointly, avoiding the independence assumptions of sequential approaches.
2. **Linear-Scale Fitting** - Operates on linear-scale PSD rather than log-scale, preserving the natural variance structure and avoiding multiplicative bias.
3. **Robust Regularization** - Elastic constraint regularization prevents overfitting and ensures physiologically plausible parameter estimates.
4. **Flexible Kernel Choices** - Supports Lorentzian (aperiodic) and Gaussian (oscillatory) kernels; can be extended with alternative kernels.
5. **Multi-Peak Detection** - Capable of detecting and characterizing multiple oscillatory peaks across frequency bands.
6. **Quality Metrics** - Provides variance-weighted R² and other diagnostics for model assessment and outlier detection.
7. **Parallel Processing** - Scales efficiently to large datasets via MATLAB's `parfor` loops.

## Citation

If you use this toolkit, please cite:

```
1. Wang, Y. et al. Whole brain resting-state EEG dynamic: A mixture of linear aperiodic and nonlinear resonant stochastic processes. Preprint at https://doi.org/10.1101/2025.06.27.661950 (2025).
2. Li, M. et al. Aperiodic and Periodic EEG Component Lifespan Trajectories: Monotonic Decrease versus Growth-then-Decline. 2025.08.26.672407 Preprint at https://doi.org/10.1101/2025.08.26.672407 (2025).
3. Li, M. et al. Harmonized-Multinational qEEG norms (HarMNqEEG). NeuroImage 256, 119190 (2022).
---

## Data link:
1-The shared raw cross-spectra with encrypted ID is hosted at
Synapse.org (10.7303/syn26712693) and complete access is possible
after login in the system. 
2-The multinational harmonized norms (HarMNqEEG norms) of traditional log-spectra and Riemannian cross-spectra
are hosted at Synapse.org (10.7303/syn26712979).
## Code link
1-The HarMNqEEG code for calculating the z-scores based on
the HarMNqEEG norm opened in GitHub, see: https://github.com/LMNonlinear/HarMNqEEG.
2-The qEEG with linearlized aperiodic and periodic norms toolbox is in https://github.com/LMNonlinear/AP-qEEG
---

## License
This project is distributed under the GNU General Public License v3.0 (see `LICENSE`).

## Authors
Ying Wang && Min Li
Oct 20, 2025




