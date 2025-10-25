function fit_bs(self)

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
self.para_bs.use_bssigma = true; %false
self.para_bs.use_bsnu = true; %false
self.para_bs.use_bsd = true; %false

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
[~, ~, ~, ~, self.para_bs] = vecpara2para(x, self.order, self.para_xi, self.para_alpha, self.para_bs, self.model, self.peak_relation);
% %}

end
