% Author: Ying Wang, Min Li
% Create Time: 2025
% Copyright (c): 2020-2025 Ying Wang, yingwangrigel@gmail.com,
%                Min Li, minli.231314@gmail.com
% Joint China-Cuba LAB, UESTC, Hangzhou Dianzi University
% License: GNU General Public License v3.0 (see LICENSE file)

function NonlinearOptimization(self)
% %{

% fit_s(self)

% new spectrum parameter
% fit_s(self);
%{
freemu = isempty(self.mufix);

if freemu
    self.mufix = self.para_alpha.kernel.para.mu;
end

fit_s(self);

if freemu
    self.mufix = [];
end
%}
% % tune spectrum para especially the center
% tune_para_s_bs(self);

fit_s(self);

% tune spectrum para especially the center
% tune_para_s_bs(self);

% bispectrum
% fit_bs(self);
% %}

% self.para_fit.FunctionTolerance=1e-6;%1e-7
% self.para_fit.OptimalityTolerance=1e-10;%1e-11 1e-6 1e-10;

% get the theoretical spectrum and bispectrum with the lost function Brillinger 1985 and Leonenko proposed in 1998
%%
%{
if ~isempty(self.bs) && size(self.bs, 1) == size(self.s, 1) && size(self.bs, 2) == size(self.s, 1)
    self.para_bs.input_config = self.para_bs; %false

    self.para_bs.use_bssigma = false; %false
    self.para_bs.use_bsnu = false; %false
    self.para_bs.use_bsd = false; %false

    if self.para_bs.fix_mu
        self.mufix = self.para_alpha.kernel.para.mu;
    end

    fit_full(self);

    % tune bispectrum
    self.para_bs.use_bssigma = self.para_bs.input_config.use_bssigma;
    self.para_bs.use_bsnu = self.para_bs.input_config.use_bsnu;
    self.para_bs.use_bsd = self.para_bs.input_config.use_bsd; %false

    if self.para_bs.use_bssigma || self.para_bs.use_bsnu || self.para_bs.use_bsd
        fit_bs(self);
    end

else
    disp('fit spectra only')
end

%}

%%{
if ~isempty(self.bs) && size(self.bs, 1) == size(self.s, 1) && size(self.bs, 2) == size(self.s, 1)
    self.para_bs.input_config = self.para_bs; %false

    self.para_bs.use_bssigma = false; %false
    self.para_bs.use_bsnu = false; %false
    self.para_bs.use_bsd = false; %false

    if self.para_bs.fix_mu
        self.mufix = self.para_alpha.kernel.para.mu;
    end

    fit_bs(self);

    % tune bispectrum
    self.para_bs.use_bssigma = self.para_bs.input_config.use_bssigma;
    self.para_bs.use_bsnu = self.para_bs.input_config.use_bsnu;
    self.para_bs.use_bsd = self.para_bs.input_config.use_bsd; %false

    if self.para_bs.use_bssigma || self.para_bs.use_bsnu || self.para_bs.use_bsd
        fit_bs(self);
    end

else
    disp('fit spectra only')
end

%%}
%%
% tic
% [self.stat.bc.xi.t, self.stat.bc.xi.ci, self.stat.bc.xi.se] = self.calc_bc_t('xi', false);
% [self.stat.bc.alpha.t, self.stat.bc.alpha.ci, self.stat.bc.alpha.se] = self.calc_bc_t('alpha', false);
%{
[self.stat.bc_polar.xi.t, self.stat.bc_polar.xi.ci, self.stat.bc_polar.xi.se] = self.calc_bc_t('xi', true);
[self.stat.bc_polar.alpha.t, self.stat.bc_polar.alpha.ci, self.stat.bc_polar.alpha.se] = self.calc_bc_t('alpha', true);
self.model.jacobian = [];
t = toc;
disp(['[calc_bc_t] time cost for t value of bicoherence estimate = ', num2str(t)]);
%}
end

