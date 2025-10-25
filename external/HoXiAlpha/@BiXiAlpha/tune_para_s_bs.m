function tune_para_s_bs(self)

% only the diagonal and cross line of fundamental peak
self.order = '1+2';

record.idx_diag = reshape((1:self.Nf) + (0:self.Nf - 1) * self.Nf, [], 1);
record.bs_diag = [self.bs(record.idx_diag)]; %self.bs(idx_mu, :)';
record.bsvar_diag = [self.bsvar(record.idx_diag)]; %self.bshatvar(idx_mu, :)';
record.bshatvar_diag = [self.bshatvar(record.idx_diag)]; %self.bshatvar(idx_mu, :)';

% x0
x0 = self.x;
x0(isnan(x0)) = 0;

% L1 = length(self.vec_para('k'))+ length(self.vec_para('xi'))+length(self.vec_para('alpha'));

% N1MAX = self.N1MAX;
% N2MAX = self.N2MAX;

% self.region_bs='mu+diag';
% lb = self.lb;
% ub = self.ub;
% % only update mu and diag of beta & gamma, fix other paras
% if strcmpi(self.para_bs.type, 'bigaussian+xi+delta')||strcmpi(self.para_bs.type, 'trigaussian+xi+delta')
%     % idx_update = reshape(find((eye(N1MAX)| [ones(N1MAX, 2),  zeros(N1MAX, N1MAX - 2)] | [zeros(N1MAX, 2), ones(N1MAX, 1), zeros(N1MAX, N1MAX - 3)] )), 1, []);
%     idx_update = reshape(find((eye(N1MAX) | [zeros(N1MAX, 2+self.ks), ones(N1MAX, 1), zeros(N1MAX, N1MAX - 3-self.ks)])), 1, []);
% elseif strcmpi(self.para_bs.type, 'bigaussian+xi') || strcmpi(self.para_bs.type, 'trigaussian+xi')
%     % idx_update = reshape(find((eye(N1MAX) | [ones(N1MAX, 1),  zeros(N1MAX, N1MAX - 1)] | [zeros(N1MAX, 1), ones(N1MAX, 1), zeros(N1MAX, N1MAX - 2)] )), 1, []);
%     idx_update = reshape(find((eye(N1MAX) | [zeros(N1MAX, 1+self.ks), ones(N1MAX, 1), zeros(N1MAX, N1MAX - 2-self.ks)])), 1, []);
% else
%     error('no such para_bs.type')
% end
%

% x0 = strict_bound_x0(x0, lb, ub);
% % idx_update = [1:L1, idx_update + L1, idx_update + N2MAX + L1,length(lb)-self.para_bs.Nsigma:length(lb)]; %, length(lb) - 2:length(lb), L1+ 2 * N2 + 1:length(x0)
% % idx_keep = setdiff(1:length(lb), idx_update);
% lb(idx_keep) = x0(idx_keep);
% ub(idx_keep) = x0(idx_keep);

self.region_bs = 'mu+diag';
lb = self.lb;
ub = self.ub;

% fit
LostObjFunTune = @(x) ObjFunBiXiAlpha(self, x, record);
% nonlcon = @(x) nonlcon_xi(x, self);
nonlcon = [];
[x] = nonlinear_fit(self, LostObjFunTune, x0, lb, ub, nonlcon);
[self.ks, self.kh, self.para_xi, self.para_alpha, self.para_bs] = vecpara2para(x, self.order, self.para_xi, self.para_alpha, self.para_bs, self.model);

%% deal with nans
self.para_bs.kernel.para.h(isnan(self.para_bs.kernel.para.h)) = 0;

% figure(11)
% plot(log10(abs([[self.s;self.bs(record.idx_mu, :)';self.bs(record.idx_diag)],[self.shat;self.bshat(record.idx_mu, :)'; self.bshat(record.idx_diag)]])))
% % figure(12)
% % surfbs(self.fx, self.fy, abs(self.bs))
% % figure(13)
% % surfbs(self.fx, self.fy, abs(self.bshat))
%

end
