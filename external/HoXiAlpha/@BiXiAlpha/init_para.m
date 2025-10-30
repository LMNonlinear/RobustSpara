% Author: Ying Wang, Min Li
% Create Time: 2025
% Copyright (c): 2020-2025 Ying Wang, yingwangrigel@gmail.com,
%                Min Li, minli.231314@gmail.com
% Joint China-Cuba LAB, UESTC, Hangzhou Dianzi University
% License: GNU General Public License v3.0 (see LICENSE file)

function init_para(self)

% self.bs = tril2full(self.bs, 0, false);
% self.bsvar = tril2full(self.bsvar, 0, false);
% self.bsdenom = tril2full(self.bsdenom, 0, false);

% if isempty(self.Rxi)
%     self.Rxi=0.5;
% end
if isempty(self.para_xi) || ~isfield(self.para_xi, 'type') || isempty(self.para_xi.type)
    self.para_xi.type = 'tstudent'; %'tstudent+peak';'tstudent'; 'linear';
end

if isempty(self.para_alpha) || ~isfield(self.para_alpha, 'type') || isempty(self.para_alpha.type)
    self.para_alpha.type = 'tstudent'; %'gaussian';
end

if isempty(self.para_bs) || ~isfield(self.para_bs, 'type') || isempty(self.para_bs.type)
    self.para_bs.type = 'tstudent'; % bigaussian+xi+delta bigaussian+xi trigaussian+xi
end

if isempty(self.para_bs) || ~isfield(self.para_bs, 'no_subharmonic_off_diag') || isempty(self.para_bs.no_subharmonic_off_diag)
    self.para_bs.no_subharmonic_off_diag = false;
end

if isempty(self.para_alpha) || ~isfield(self.para_alpha, 'no_zero_peak') || isempty(self.para_alpha.no_zero_peak)
    % self.para_alpha.no_zero_peak = true;
    self.para_alpha.no_zero_peak = false;
end

if isempty(self.para_bs) || ~isfield(self.para_bs, 'no_xi_alpha') || isempty(self.para_bs.no_xi_alpha)
    % warning('no_xi_alpha is set to false')
    % self.para_bs.no_xi_alpha = false;
    self.para_bs.no_xi_alpha = true;
end

if isempty(self.para_bs) || ~isfield(self.para_bs, 'fix_fxy') || isempty(self.para_bs.fix_fxy)
    self.para_bs.fix_fxy = true;
end

if isempty(self.para_xi) || ~isfield(self.para_xi, 'use_bsd') || isempty(self.para_xi.use_bsd)
    self.para_xi.use_xid = true; %false
    % self.para_xi.use_xid = false; %false
end

if isempty(self.para_bs) || ~isfield(self.para_bs, 'fix_mu') || isempty(self.para_bs.fix_mu)
    % self.para_bs.use_bssigma = false; %false
    self.para_bs.fix_mu = true; %false
end

if isempty(self.para_bs) || ~isfield(self.para_bs, 'use_bssigma') || isempty(self.para_bs.use_bssigma)
    % self.para_bs.use_bssigma = false; %false
    self.para_bs.use_bssigma = true; %false
end

if isempty(self.para_bs) || ~isfield(self.para_bs, 'use_bsnu') || isempty(self.para_bs.use_bsnu)
    % self.para_bs.use_bsnu = false; %false
    self.para_bs.use_bsnu = true; %false
end

if isempty(self.para_bs) || ~isfield(self.para_bs, 'use_bsd') || isempty(self.para_bs.use_bsd)
    % self.para_bs.use_bsd = false; %false
    self.para_bs.use_bsd = true; %false
end

if isempty(self.para_bs) || ~isfield(self.para_bs, 'use_bsb') || isempty(self.para_bs.use_bsb)
    self.para_bs.use_bsb = false; %false
    % self.para_bs.use_bsb = true; %false
end

if isempty(self.para_bs) || ~isfield(self.para_bs, 'additive') || isempty(self.para_bs.additive)
    self.para_bs.additive = true;
end

if isempty(self.para_alpha) || ~isfield(self.para_alpha, 'mu_df_lb') || isempty(self.para_alpha.mu_df_lb)
    self.para_alpha.mu_df_lb = 1;
end

if isempty(self.para_alpha) || ~isfield(self.para_alpha, 'mu_df_ub') || isempty(self.para_alpha.mu_df_ub)
    self.para_alpha.mu_df_ub = 1;
end

if isempty(self.para_fit) || ~isfield(self.para_alpha, 'optfun') || isempty(self.para_fit.optfun)
    self.para_fit.optfun = 'lsqnonlin'; %fmincon %lsqnonlin

    if strcmpi(self.para_fit.optfun, 'fmincon')
        self.para_fit.Algorithm = 'interior-point'; % 'interior-point' 'levenberg-marquardt' 'trust-region-reflective' 'sqp' 'active-set'
    else
        self.para_fit.Algorithm = 'levenberg-marquardt'; %;'levenberg-marquardt'; % 'trust-region-reflective', 'levenberg-marquardt', or 'interior-point'
    end

    % self.para_fit.MaxIterations = 10000;%10000;%120000;%6000
    % self.para_fit.MaxFunctionEvaluations =80000;%40000;
    % self.para_fit.Display = 'none';%'final-detailed'; %'final-detailed';%'none',  'off',  'iter',  'iter-detailed', 'notify',  'notify-detailed',  'final',  or 'final-detailed'
    % self.para_fit.UseParallel=false;%~isinparfor();%~isinparfor();%~isinparfor();%warning('close parfor');~isinparfor();
    % self.para_fit.FunctionTolerance=1e-8;
    % self.para_fit.OptimalityTolerance=1e-12;%1e-6 1e-10;

    if self.Ns > 1
        self.para_fit.MaxIterations = 500; %1000 3000; %1500; %10000; %120000; %6000 4000 8000
        self.para_fit.MaxFunctionEvaluations = 2000; %4000 12000; %6000; %40000; 16000 32000
    else
        self.para_fit.MaxIterations = 2000; %1000 3000; %1500; %10000; %120000; %6000 4000 8000
        self.para_fit.MaxFunctionEvaluations = 8000; %4000 12000; %6000; %40000; 16000 32000
    end

    self.para_fit.Display = 'none'; %'final-detailed'; %'final-detailed'; %'none','off', 'iter', 'iter-detailed', 'notify', 'notify-detailed', 'final', or 'final-detailed'
    self.para_fit.UseParallel = false; %~isinparfor(); false %~isinparfor(); %warning('close parfor'); ~isinparfor();
    self.para_fit.ScaleProblem = 'jacobian'; %'none
    % if strcmpi(self.para_fit.Algorithm, 'levenberg-marquardt')
    % self.para_fit.ScaleProblem = 'jacobian'; %'Jacobian'; %'none
    % else
    % self.para_fit.ScaleProblem = true;
    % end
    self.para_fit = set_defaults(self.para_fit, 'StepTolerance', 1e-8);
    self.para_fit = set_defaults(self.para_fit, 'FunctionTolerance', 1e-8);
    self.para_fit = set_defaults(self.para_fit, 'OptimalityTolerance', 1e-8);
    self.para_fit = set_defaults(self.para_fit, 'FiniteDifferenceStepSize', 1e-8);

    % self.para_fit.StepTolerance = 1e-8; %1e-8; 1e-8; %1e-6;
    % self.para_fit.FunctionTolerance = 1e-8; %1e-8; 1e-8; %1e-6; %1e-7 1e-8
    % self.para_fit.OptimalityTolerance = 1e-8; %1e-8; 1e-12; %1e-10; %1e-11 1e-6 1e-10; 1e-12
    % self.para_fit.FiniteDifferenceStepSize = 1e-8; %1e-6

    % self.para_fit.MaxIterations = 1000;%10000;%120000;%6000
    % self.para_fit.MaxFunctionEvaluations =4000;%40000;
    % self.para_fit.Display = 'none';%'final-detailed'; %'final-detailed';%'none',  'off',  'iter',  'iter-detailed', 'notify',  'notify-detailed',  'final',  or 'final-detailed'
    % self.para_fit.UseParallel=false;%~isinparfor();%~isinparfor();%~isinparfor();%warning('close parfor');~isinparfor();
    % self.para_fit.FunctionTolerance=1e-5;
    % self.para_fit.OptimalityTolerance=1e-9;%1e-6 1e-10;
end

if ~isempty(self.s)

    if isempty(self.svar)
        self.svar = ones(size(self.s, 1), 1);
    end

    % xi proccess
    init_xi(self);
    % alphas proccess
    init_alpha(self);

end

if ~isempty(self.bs)

    % idxused for the square matrix
    get_bifrequency(self);
    self.tfused_segs = ~isnan(self.bsdenom) & ~isinf(self.bsdenom) & self.bsdenom ~= 0 & ...
        ~isnan(self.bs) & ~isinf(self.bs) & self.bs ~= 0;
    self.tfused_segs_sep_z = repmat(self.tfused_segs, [1, 1, 2]);

    self.tfused = logical(prod(self.tfused_segs, 3));
    self.idxused = find(logical(self.tfused));
    % bispectrum
    init_bs(self);
    init_p_q(self);

end

end

function init_p_q(self)

% %{
r_s = mean(abs(calc_loss_s(self, self.loss_type, self.s, mean(self.s, 1, "omitnan"), self.svar)), "all", "omitnan");
% bshatvar = sqrt(get_sxsysxy(self.f, self.s, [], self.fxy));
bsvar = self.bshatdenom;
r_bs = mean(abs(calc_loss_bs(self, self.loss_type, self.bc, mean(self.bs, [1, 2], "omitnan"), [], bsvar)), "all", "omitnan");
% r_bs = mean(abs(calc_loss_bs(self, self.loss_type, self.bc, mean(self.bs, [1, 2],"omitnan"), [], self.bshatvar)),"all","omitnan");
% r_s = sqrt(mean((self.s - mean(self.s, 1,"omitnan")) .^ 2,"omitnan"));
% diff_bs=self.bs - mean(self.bs, [1, 2],"omitnan");
% diff_bs=[real(diff_bs);imag(diff_bs)];
% r_bs = sqrt(mean((diff_bs) .^ 2, [1, 2],"omitnan"));
self.p = r_bs ./ (r_s + r_bs);
self.q = r_s ./ (r_s + r_bs);
% %}

r_s = mean(self.svar, 1);
r_bs = mean(self.bsvar(:), 1, "omitnan");
self.p = r_bs ./ (r_s + r_bs);
self.q = r_s ./ (r_s + r_bs);

end

