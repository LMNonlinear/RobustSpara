% Author: Ying Wang, Min Li
% Create Time: 2025
% Copyright (c): 2020-2025 Ying Wang, yingwangrigel@gmail.com,
%                Min Li, minli.231314@gmail.com
% Joint China-Cuba LAB, UESTC, Hangzhou Dianzi University
% License: GNU General Public License v3.0 (see LICENSE file)

function fit_alpha(self)
% increase accuracy of this fit
% warning('increase accuracy of this fit');

% para_fit=self.para_fit;
% self.para_fit.MaxIterations = 3000;%1500;%10000;%120000;%6000 4000 8000
% self.para_fit.MaxFunctionEvaluations =12000;%6000;%40000; 16000 32000
% self.para_fit.StepTolerance=1e-8;%1e-6;
% self.para_fit.FunctionTolerance=1e-8;%1e-6;%1e-7 1e-8
% self.para_fit.OptimalityTolerance=1e-12;%1e-10;%1e-11 1e-6 1e-10; 1e-12

% fmincon
% spectrum
self.order = '1';
self.iiter = 0;
lb = self.lb;
ub = self.ub;
x0 = self.x;

lb((length(bound_k(self)) + 1):(length(bound_k(self)) + 1) + length(bound_xi(self))) = x0((length(bound_k(self)) + 1):(length(bound_k(self)) + 1) + length(bound_xi(self)));
ub((length(bound_k(self)) + 1):(length(bound_k(self)) + 1) + length(bound_xi(self))) = x0((length(bound_k(self)) + 1):(length(bound_k(self)) + 1) + length(bound_xi(self)));

% nonlcon = @(x) nonlcon_xi(x, self);
nonlcon = [];
[x] = nonlinear_fit(self, @(x) ObjFunBiXiAlpha(self, x), x0, lb, ub, nonlcon);
[self.ks, self.kh, ~, self.para_alpha] = vecpara2para(x, self.order, self.para_xi, self.para_alpha, self.para_bs, self.model);

end

