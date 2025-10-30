% Author: Ying Wang, Min Li
% Create Time: 2025
% Copyright (c): 2020-2025 Ying Wang, yingwangrigel@gmail.com,
%                Min Li, minli.231314@gmail.com
% Joint China-Cuba LAB, UESTC, Hangzhou Dianzi University
% License: GNU General Public License v3.0 (see LICENSE file)

function fit_xi(self)

x0 = self.para_xi.kernel.vpara;

h = 0.4;
[s] = LLSmooth(self.f, mean(self.s, 2), h, self.f);

% s = exp(mean(log(self.s), 2));

[lb, ub] = bound_xi(self);
% lb(1) = 1/10 .* max(self.trend);
% ub(1) = 10 .* max(self.trend);

x = nonlinear_fit(self, @(x) LostObjFun(x, self), x0, lb, ub);
self.para_xi.kernel.vpara = x;
self.trend=pred_s_xi(self.f, self.para_xi);
end
%% lost function for the xi process
function loss = LostObjFun(x, self)
    x = recover_fixed_para(x, [], [], self.para_fit.parainfo);
    para_xi = self.para_xi;

    para_xi.kernel.vpara = x;

    shat_xi = pred_s_xi(self.f, para_xi);
    self.iiter = self.iiter + 1;
    % plot(abs([self.trend,xihat]))
    if self.verbose && mod(self.iiter, 100) == 1
        figure(31)
        plot(self.f, log10(abs([self.s, shat_xi])))
        drawnow
    end

    loss = calc_loss_s(self, self.loss_type, self.s, shat_xi, self.svar, []);
    diff_endpoint = [self.s(1) - shat_xi(1); self.s(end) - shat_xi(end)];
    trend_residual = shat_xi - self.s;
    trend_residual(trend_residual < 0) = 0;
    
    penalty = [ 1e3.* diff_endpoint(1); ... %self.lambda(3)
                2*1e3.*diff_endpoint(2);...
                   1e4 .* trend_residual];%self.lambda(4)
    loss = [loss; penalty];
    self.model.loss_history(self.iiter, :) = sum(loss(~isnan(loss)));
end

%{
idx = (self.f < self.para_xi.f_max);
self.trend(idx) = s(idx);
% also need avoid alpha peak in this range
[~, ~, trend, ~, stats] = RobustDetrend(log(self.trend), 2, 0.975, self.f);
idx_peak = (abs(stats.resid) > stats.mad_s & stats.resid > 0);
self.trend(idx_peak) = exp(trend(idx_peak));

% remove the upflip of whole spectrum, due to the robust fit
idx = (self.f > self.para_xi.f_max);
dif = (log(self.trend) - log(s));
% mad_dif = median(abs(dif - median(dif)));
shift = max(dif(dif > 0)); % & abs(dif)<mad_dif
self.trend(idx) = exp(log(self.trend(idx)) - shift);
% plot(log10([self.trend]))

[lb, ub] = bound_xi(self);
% lb(1) = 1/10 .* max(self.trend);
% ub(1) = 10 .* max(self.trend);

[x] = nonlinear_fit(self, @(x) LostObjFun(x, self), x0, lb, ub);

%% save result
self.para_xi.kernel.vpara = x;

end

%% lost function for the xi process
function loss = LostObjFun(x, self)
x = recover_fixed_para(x, [], [], self.para_fit.parainfo);
para_xi = self.para_xi;

para_xi.kernel.vpara = x;

xihat = pred_s_xi(self.f, para_xi);
self.iiter = self.iiter + 1;
% plot(abs([self.trend,xihat]))
if self.verbose && mod(self.iiter, 20) == 1
    figure(31)
    plot(self.f, log10(abs([self.trend, xihat])))
    drawnow
end

% loss = sum(log(xihat)+self.trend./xihat);% loss=(xihat-trend).^2;
% loss = calc_loss_s(self, self.loss_type, self.trend, xihat, ones(size(xihat)), self.svar);
% loss = calc_loss_s(self, self.loss_type, self.trend, xihat, ones(size(xihat)));
loss = calc_loss_s(self, self.loss_type, self.trend, xihat, self.svar, []);
self.model.loss_history(self.iiter, :) = sum(loss(~isnan(loss)));
end
%}

