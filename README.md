# RobustXiAlpha

RobustXiAlpha is a MATLAB toolbox for robustly decomposing quantitative EEG (qEEG) power spectra into aperiodic (``\xi``) and oscillatory (alpha-band) constituents. It wraps the `BiXiAlpha` class and several companion utilities that smooth spectra, initialize non-linear fits, and export parameter tables for downstream normative modelling.

## Highlights
- Implements Lorentzian aperiodic and Gaussian oscillatory kernels that are jointly optimized with Levenberg-Marquardt for stable estimates of the ``\xi`` slope, offsets, peak frequency, bandwidth, and power.
- Provides scripts that process HarMNqEEG cross-spectra and export channel-wise tables of alpha and ``\xi`` parameters for normative analyses.
- Ships quality metrics (e.g., weighted R^2) and helpers for smoothing, bounds handling, bispectral statistics, and table manipulation.
- Ready to integrate with broader qMEEG pipelines by reusing the existing QMEEG utilities (e.g., `get_subtable`, `test_folder`) and parallel execution helpers (`startmatlabpool`).

## Methodological Context
- **Multinational data workflow** - As described in *xialpha_Age_V4_method_biorxiv2.docx*, HarMNqEEG recordings originate from 14 studies conducted in 9 countries on 12 EEG devices. Data contributors run a shared script locally, upload only cross-spectra plus metadata (age, sex, country, device), and the central site harmonizes batches after anonymisation via Arng2. RobustXiAlpha's examples mirror this decentralised ingest.
- **Robust spectral parametrisation** - Following the same method note, the aperiodic background is modelled with a Lorentzian kernel (height ``h``, exponent ``\chi``, knee, and baseline), while periodic peaks are Gaussian (centre ``\mu``, scale ``\sigma``, amplitude ``h``). Both components are fitted simultaneously to the linear-scale PSD, avoiding the multiplicative assumptions made by log-domain approaches such as FOOOF.

## Repository Layout
- `utility/` - Core helpers: smoothing kernels, bounds classification, Jacobian estimation, conversions, plotting, table utilities, etc.
- `external/HoXiAlpha/` - The `BiXiAlpha` class definition and optimisation routines shared with the HoXiAlpha project.
- `external/` *(others)* - Third-party resources used by the toolbox (parallel helpers, colormaps, detrending, brain eigenmodes).
- `example/qMEEG/` - End-to-end scripts demonstrating how to read HarMNqEEG cross-spectra, run RobustXiAlpha fits, and export summary tables.

## Requirements
- MATLAB R2021a or newer (tested with >= R2022b).
- Signal Processing Toolbox (for PSD utilities) and Parallel Computing Toolbox (scripts use `parfor`).
- Access to HarMNqEEG spectral tables or your own PSD matrices with compatible metadata.
- The broader QMEEG utility set on the MATLAB path (for functions such as `test_folder`).

## Installation & Setup
1. Clone or download `RobustXiAlpha` into your MATLAB workspace.
2. Add the toolbox folders to the MATLAB path:
   ```matlab
   addpath(genpath(fullfile(pwd, 'utility')));
   addpath(genpath(fullfile(pwd, 'external')));
   addpath(genpath(fullfile(pwd, 'example')));
   ```
3. Ensure the QMEEG shared utilities that provide `test_folder`, `get_folders_in_directory`, etc., are also on the path.
4. If you plan to use parallel execution, confirm that a MATLAB parallel pool can be opened (`parpool`, or the provided `startmatlabpool` wrapper).

## Quick Start
1. **Prepare cross-spectra** - Follow the HarMNqEEG workflow or adapt your own dataset into tables matching `example/qMEEG/data/HarMNqEEG/ystarlog_1563.csv`, with frequency (`freq`) columns plus channel-wise log power entries.
2. **Run the smoothing + fit demo** - Execute `demo_xialpha_ystar_smooth_lorenz_gaussian.m`. The script:
   - Boots a parallel pool (up to 20 workers).
   - Reads a subject's cross-spectrum, smooths it via local Gaussian kernels (`NwSmoothInline`), and assembles the `hos` structure (`f`, `s`, `Fs`, metadata).
   - Instantiates `BiXiAlpha` objects for each EEG channel, switches on robust elastic constraints, and fits both aperiodic and alpha components.
   - Saves channel-wise results under `result/qMEEG/<task_name>/subj_<id>.mat` and renders diagnostic plots (`xialpha(ichan).show`).
3. **Extract developmental trajectories** - Run `demo_xialpha_ystar_trajectory_1.m` on the saved `.mat` files. It iterates over all channels, collates per-band alpha peaks (delta/theta/alpha/beta), exports:
   - `T_alpha_hoxialpha*.csv` with every detected peak,
   - `T_xi_hoxialpha*.csv` with the Lorentzian exponents/offsets,
   - `T_alpha_band_hoxialpha*.csv` with the strongest peak per canonical band.

## Input Specification
The main entry point (`BiXiAlpha`) expects:
- `f` - frequency vector (Hz) sampled uniformly, typically <= 20 Hz for HarMNqEEG.
- `s` - power spectral density values (linear scale), columns per channel.
- `bs` *(optional)* - bispectral terms; pass `[]` if not available.
- `ks`, `kh` - taper parameters for spectral smoothing (default 4 / 8 in demos).
- `Fs` - sampling rate (Hz) used for normalisation.
Metadata (age, sex, country, device, dataset) are attached to each `BiXiAlpha` instance via the `.info` struct.

## Outputs
For every channel, `BiXiAlpha` stores:
- `para_xi.kernel.para` - Lorentzian parameters (`chi`, `h`, `mu`, offsets) governing the aperiodic ``\xi`` component.
- `para_alpha.kernel.para` - Gaussian peak parameters (`mu`, `sigma`, `h`) for the dominant alpha rhythm.
- `r2`, `r2w`, `r2_s` - Quality-of-fit metrics (plain and variance-weighted).
- `loss`, `regularization`, `lambda` - Optimisation settings used to reach the solution.
Saved `.mat` files preserve arrays of `BiXiAlpha` objects that can be reloaded for additional analysis, plotting, or band-specific summarisation.

## Quality Control & Harmonisation
- Use the weighted R^2 (`calc_r2w`) to compare model and empirical spectra when observation variances are known.
- The method note recommends screening outliers using global z-scores derived from age-adjusted normative surfaces and harmonising batches with ComBat-like adjustments before fitting; apply those steps upstream when processing multi-site data.

## Extending the Toolbox
- Replace or augment `psd_smooth_nw` with alternative smoothers (`LLSmooth`, Savitzky-Golay) as needed.
- Toggle `xialpha(:,ich).peak_relation`, `regularization`, or `lambda` to enforce anatomical priors (e.g., harmonic constraints) or relax elastic penalties.
- Batch-processing scripts can be adapted to other montages by updating channel headers (`get_spec_hearder`) and metadata handling.

## License
This project is distributed under the GNU General Public License v3.0 (see `LICENSE`).

## Authors
Ying Wang && Min Li

