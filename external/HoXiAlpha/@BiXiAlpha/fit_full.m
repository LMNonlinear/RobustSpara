function fit_full(self)

% disp(['start fit for name: [',num2str(self.name),']']);

% avoid redundant
record.shat = self.shat;
record.bs = self.bs;
record.bc = self.bc;

switch self.loss_type
    case {'bicoherence'}
        record.bsvar = self.bsdenom;
        record.bshatvar = self.bshatdenom;
    case {'HeteroEstVar', 'HeteroscedasticityEstVar'}
        record.bsvar = self.bsvar;
        record.bshatvar = self.bshatvar;
end

% %{
% warning('use init val')
% bispectrum
self.order = '2';
self.region_bs = 'all';
lb = self.lb;
ub = self.ub;
x0 = self.x;
% x0 = strict_bound_x0(x0, lb, ub);
% idx_update = 1:length(lb);
% idx_keep = setdiff(1:length(lb), idx_update);
% lb(idx_keep) = x0(idx_keep);
% ub(idx_keep) = x0(idx_keep);

[x] = nonlinear_fit(self, @(x) ObjFunBiXiAlpha(self, x, record), x0, lb, ub);
[~, ~, ~, ~, self.para_bs] = vecpara2para(x, self.order, self.para_xi, self.para_alpha, self.para_bs, self.model,self.peak_relation);
% %}

% spectrum and bispectrum
% %{
self.order = '1+2';
self.region_bs = 'all';
lb = self.lb;
ub = self.ub;
x0 = self.x;

[x] = nonlinear_fit(self, @(x) ObjFunBiXiAlpha(self, x, record), x0, lb, ub);
[self.ks, self.kh, self.para_xi, self.para_alpha, self.para_bs] = vecpara2para(x, self.order, self.para_xi, self.para_alpha, self.para_bs, self.model,self.peak_relation);
% self.loss_type=loss_type;
% %}

% %%
% tic
% % [self.stat.bc.xi.t, self.stat.bc.xi.ci, self.stat.bc.xi.se] = self.calc_bc_t('xi', false);
% % [self.stat.bc.alpha.t, self.stat.bc.alpha.ci, self.stat.bc.alpha.se] = self.calc_bc_t('alpha', false);

% [self.stat.bc_polar.xi.t, self.stat.bc_polar.xi.ci, self.stat.bc_polar.xi.se] = self.calc_bc_t('xi', true);
% [self.stat.bc_polar.alpha.t, self.stat.bc_polar.alpha.ci, self.stat.bc_polar.alpha.se] = self.calc_bc_t('alpha', true);
% self.model.jacobian = [];
% t = toc;
% disp(['[calc_bc_t] time cost for t value of bicoherence estimate = ', num2str(t)]);

self.model.jacobian = [];

end
