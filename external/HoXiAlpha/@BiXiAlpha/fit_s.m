% Author: Ying Wang, Min Li
% Create Time: 2025
% Copyright (c): 2020-2025 Ying Wang, yingwangrigel@gmail.com,
%                Min Li, minli.231314@gmail.com
% Joint China-Cuba LAB, UESTC, Hangzhou Dianzi University
% License: GNU General Public License v3.0 (see LICENSE file)

function fit_s(self)
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


% nonlcon = @(x) nonlcon_xi(x, self);
nonlcon = [];
[x] = nonlinear_fit(self, @(x) ObjFunBiXiAlpha(self, x), x0, lb, ub, nonlcon);
[self.ks, self.kh, self.para_xi, self.para_alpha] = vecpara2para(x, self.order, self.para_xi, self.para_alpha, self.para_bs, self.model,self.peak_relation);

% self.para_fit.MaxIterations=para_fit.MaxIterations;
% self.para_fit.MaxFunctionEvaluations=para_fit.MaxFunctionEvaluations;
% self.para_fit.StepTolerance=para_fit.StepTolerance;
% self.para_fit.FunctionTolerance=para_fit.FunctionTolerance;
% self.para_fit.OptimalityTolerance=para_fit.OptimalityTolerance;

end

