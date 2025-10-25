%
% Main script for robust spectral parameter estimation
%
% This script demonstrates how to use the BiXiAlpha class to decompose
% quantitative EEG (qEEG) power spectra into aperiodic (xi) and oscillatory
% (alpha-band) components.
%
% Input:
%   - Spec data (47x1 column vector) from data/Spec.mat
%   - Log-scale power spectral density values
%
% Output:
%   - BiXiAlpha model results saved in results/robust_s_para/robust_s_para_analysis/
%   - xialpha_result.mat containing fitted BiXiAlpha objects
%
% Usage:
%   1. Set debug=true to process only the first channel (for testing)
%   2. Set debug=false to process all channels (requires parallel computing)
%   3. Run the script: run('example/main_robust_s_para.m')
%   4. Results will be saved in the results folder
%
% Key parameters:
%   - ks: Number of harmonics for xi (aperiodic) component (default: 4)
%   - kh: Number of harmonics for alpha (oscillatory) component (default: 8)
%   - regularization: Type of regularization (default: elastic_constraint_mu_sigma_chi_xi)
%   - lambda: Regularization weights for different parameters
%

clear; clc;
close all

% Change to project root directory
cd(fileparts(fileparts(mfilename('fullpath'))));

% Setup paths - add all necessary directories
restoredefaultpath;
addpath(genpath('./utility/'));
addpath(genpath('./external/'));
addpath('./external/nnda/utility/');

% Setup parallel computing
maxworker = 20;
startmatlabpool(min(feature('numcores'), maxworker));

% Debug mode: set to true for single channel processing, false for all channels
debug = true;

% Setup paths
path_save = './results/robust_s_para/';
task_name = 'robust_s_para_analysis';

% Load spectral data
Spec=readtable('./data/Spec_example.csv');
% Prepare data structure
% Spec is a 47x1 column vector representing log power spectral density
% We need to convert it to linear scale and set up frequency vector and sampling rate
hos.s = exp(Spec.ystarlog9_9);  % Convert from log scale to linear scale
hos.f = Spec.freq;  % Frequency vector (1-47 Hz)
hos.Fs = max(hos.f) * 2;  % Sampling frequency (2 * max frequency)
hos.debug = debug;

% Apply smoothing to the spectrum
hos = psd_smooth_nw(hos);

% Extract frequency and spectrum
[f, s, Fs] = deal(hos.f, hos.s, hos.Fs);

% BiXiAlpha fitting parameters
ks = 4;  % Number of harmonics for xi component
kh = 8;  % Number of harmonics for alpha component

% Determine which channels to process
if hos.debug
    id_ch = 1;  % Process only first channel in debug mode
else
    id_ch = 1:size(s, 2);  % Process all channels
end

% Initialize BiXiAlpha models
for ich = id_ch
    xialpha(ich) = BiXiAlpha(f, s(:, ich), [], ks, kh, Fs);
    xialpha(ich).verbose = hos.debug;
    xialpha(ich).multiple_seg = false;

    % Set xi component parameters (Lorentzian)
    xialpha(ich).para_xi.type = 'lorentzian';

    % Set alpha component parameters (Gaussian)
    xialpha(ich).para_alpha.no_zero_peak = true;
    xialpha(ich).para_alpha.type = 'gaussian';

    % Set regularization and lambda parameters
    xialpha(ich).regularization = 'elastic_constraint_mu_sigma_chi_xi';
    xialpha(ich).lambda = [1; 1; 1e1; 1e1; 1e3];  % mu sigma startend noneg

    % Set peak relation and fitting options
    xialpha(ich).peak_relation = 'free';
    xialpha(ich).para_fit.StepTolerance = 1e-12;
    xialpha(ich).para_fit.FunctionTolerance = 1e-12;
    xialpha(ich).para_fit.OptimalityTolerance = 1e-12;
    xialpha(ich).para_fit.FiniteDifferenceStepSize = 1e-8;
end

% Fit the models
if hos.debug
    for ich = id_ch
        tempmodel = xialpha(ich);
        tempmodel.fit;
        tempmodel.temp = [];
        tempmodel.model = [];
        xialpha(ich) = tempmodel;
    end
else
    parfor ich = id_ch
        tempmodel = xialpha(ich);
        tempmodel.fit;
        tempmodel.temp = [];
        tempmodel.model = [];
        xialpha(ich) = tempmodel;
    end
end

% Save results
file_save = fullfile(path_save, task_name, 'xialpha_result.mat');
test_folder(file_save);
save(file_save, 'xialpha', '-v7.3');

% Display results
if hos.debug
    xialpha(id_ch).show;
end

disp('Analysis completed successfully!');
disp(['Results saved to: ', file_save]);

%% Helper functions

function hos = psd_smooth_nw(hos)
% Smooth the power spectral density using Nadaraya-Watson kernel
f_interp = linspace(min(hos.f), max(hos.f), length(hos.f) * 3)';
h = 0.3;  % Bandwidth parameter

hos.s = permute(hos.s, [1, 3, 2]);
[hos.s] = NwSmoothInline(hos.f, hos.s, h, f_interp);
hos.f = f_interp;
hos.s = permute(hos.s, [1, 3, 2]);
end

function [yq, L, dbg] = NwSmoothInline(x, y, h, xq)
% Nadaraya-Watson kernel smoothing

if nargin < 4 || isempty(xq)
    xq = x;
end

[n, dx] = size(x);
[nq, ~] = size(xq);

if size(x, 2) > 1
    x = reshape(x, [n, 1, dx]);
end

if size(xq, 2) > 1
    xq = reshape(xq, [nq, 1, dx]);
end

D = x - permute(xq, [2, 1, 3]);
Dkn = gaussian_kernel(D, h);
yq = permute(sum(Dkn .* y, 1) ./ sum(Dkn, 1), [2, 1, 3]);
dbg.s = sum(Dkn, 1);

if nargout > 1
    L = Dkn ./ sum(Dkn, 1).';
end
end

function u = gaussian_kernel(u, b)
% Gaussian kernel function
u = u ./ b;
u = (1 / sqrt(2 * pi)) * exp(-0.5 * sum((u .^ 2), 3));
end

function test_folder(file_path)
% Create folder if it doesn't exist
[folder, ~, ~] = fileparts(file_path);
if ~isdir(folder)
    mkdir(folder);
end
end

