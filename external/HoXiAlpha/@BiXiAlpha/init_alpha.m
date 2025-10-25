function init_alpha(self)

%% get init value
init_para_alpha(self);
bound_alpha(self);
range = [self.para_alpha.kernel.lb.mu, self.para_alpha.kernel.ub.mu];

if strcmpi(self.method_alpha, 'findmax') || (strcmpi(self.peak_relation, 'free')) %isempty(self.bs) &&
    h = 0.03;
    [s] = NwSmooth(self.f, mean(self.s, 2), h, self.f);
    % [s] = LLSmooth(self.f, mean(self.s, 2), h, self.f);

    min_distance = 0.5; % Minimum distance between peaks in frequency units
    % [mu, h, baseline, s_detrended] = find_largest_k_peaks(self.f, s, self.k, min(self.f), max(self.f) - self.para_alpha.mu_df_ub, min_distance);

    [mu, h, baseline, s_detrended] = find_largest_k_peaks(self.f, s, self.k, min(self.f) + self.para_alpha.mu_df_lb, max(self.f) - self.para_alpha.mu_df_ub, min_distance);
    h(h < 0) = 0;

    % Find the index of the peak with the highest height
    [~, idx_max_h] = max(h);

    % Frequency of the peak with the highest height
    f_max_h = mu(idx_max_h);
    if ~self.para_alpha.no_zero_peak
        mu=[0;mu];
        h=[0;h];
        % Count the number of peaks with frequency smaller than f_max_h
        self.ks = sum(mu < f_max_h)-1;
        self.ksmax = self.ks;
    else
        % Count the number of peaks with frequency smaller than f_max_h
        self.ks = sum(mu < f_max_h);
        self.ksmax = self.ks;
    end

    % Count the number of peaks with frequency larger than f_max_h
    self.kh = sum(mu >= f_max_h);
    self.khmax = self.kh;
    % k_found=length(self.para_alpha.kernel.para.mu);
    % self.k=k_found;
    init_para_alpha(self);
    bound_alpha(self);
    self.para_alpha.kernel.para.mu = mu;
    % self.para_alpha.kernel.para.h=h;
    shat_xi = pred_s_xi(self.f, self.para_xi);
    % self.para_alpha.kernel.para.h=self.para_alpha.kernel.para.h+shat_xi(ismember(self.f,mu));
    s_alpha = self.s - shat_xi;
    s_alpha_min = min(s_alpha);
    s_alpha = s_alpha - s_alpha_min;
    self.para_alpha.kernel.para.h = s_alpha(ismember(self.f, mu));

    if self.verbose
        figure(191); clf; subplot(2, 1, 1); hold on; plot(self.f, [self.s, shat_xi, s_alpha]);
        plot(mu, [self.para_alpha.kernel.para.h], 'ro', 'MarkerSize', 8, 'LineWidth', 2);
        legend(['s', 'shat_xi', 's_alpha']);
        figure(191); subplot(2, 1, 2); plot(self.f, [self.s, s, baseline, s_detrended]); hold on;
        baseline_at_mu = interp1(self.f, baseline, self.para_alpha.kernel.para.mu);
        plot(self.para_alpha.kernel.para.mu, h + baseline_at_mu, 'ro', 'MarkerSize', 8, 'LineWidth', 2);
        legend(['s', 's_smooth', 'baseline', 's-detrended']);
        % figure(191);clf;plot(self.f,10*log10([self.s,s,baseline,s_detrended]+1e-10));hold on;
        % baseline_at_mu = interp1(self.f, baseline, self.para_alpha.kernel.para.mu);
        % plot(self.para_alpha.kernel.para.mu, 10*log10(self.para_alpha.kernel.para.h + baseline_at_mu+1e-10), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
    end

    self.para_alpha.mu_init = mu;
    self.para_alpha.sigma_init = 0.8;
    self.fx = [];
    self.fy = [];
    self.fxy = [];
    self.fxfy = [];
    self.para_bs.idx_fxfy = [];
    % self.sxy=[];
    % self.sxyhat=[];
    % self.sxsysxy=[];
    % self.sxsysxyhat=[];
    % self.bsdenom=[];
    % self.bshatdenom=[];
    % self.bshatvar=[];
elseif strcmpi(self.method_alpha, 'findmax') || strcmpi(self.peak_relation, 'harmonic')
    h = 0.6;
    [s] = NwSmooth(self.f, mean(self.s, 2), h, self.f);

    [~, self.para_alpha.kernel.para.mu] = findmax_spectrum(log10(s - pred_s_xi(self.f, self.para_xi)), self.f, range); %
    s_alpha = exp(mean(log(s), 2)) - pred_s_xi(self.f, self.para_xi);
    s_alpha_min = min(s_alpha);
    s_alpha = s_alpha - s_alpha_min;
    opt.maxopt = 'raw';
    % self.para_alpha.kernel.para.h = find_max_around_xt(self.f, s_alpha, self.para_alpha.kernel.para.mu * (1:self.kh)', min(self.width_candidate, self.para_alpha.kernel.para.mu / 2), opt) + s_alpha_min;
    self.para_alpha.kernel.para.h = find_max_around_xt(self.f, s_alpha, get_harmonic_mu_sigma(self.para_alpha.kernel.para.mu, 1, self.ks, self.kh, self.para_alpha.no_zero_peak)', min(self.width_candidate, self.para_alpha.kernel.para.mu / 2), opt) + s_alpha_min;
    self.para_alpha.kernel.para.h = self.para_alpha.kernel.para.h(:);
    self.para_alpha.kernel.para.h(1:self.k0 + self.ks) = 1e-1 * self.para_alpha.kernel.para.h(1:self.k0 + self.ks);
    self.para_alpha.isubharmonic = false;
elseif strcmpi(self.method_alpha, 'maxbs')
    % tempbcnormalization = self.normalization;
    % self.normalization = 'skewness';

    if isempty(self.mufix)
        h = 0.6;
        [s] = NwSmooth(self.f, mean(self.s, 2), h, self.f)

        h = 0.6;
        bs_diag = abs(diag(mean(self.bs, 3)));
        bs_diag(isinf(abs(bs_diag))) = nan;
        bs_diag = NwSmooth(self.f, bs_diag, h, self.f)

        h = 0.6;
        bc_diag = abs(diag(mean(self.bs, 3)));
        bc_diag(isinf(abs(bc_diag))) = nan;
        bc_diag = NwSmooth(self.f, bc_diag, h, self.f)

        idx = find(self.f > range(1));
        bs_diag = bs_diag(idx);
        [maxbs, idx_max] = max(bs_diag);

        if idx_max == 1
            % might be the chaotic signal
            bc_diag = bc_diag(idx);
            [maxbc, idx_max] = max(bc_diag);
        end

        self.para_alpha.kernel.para.mu = self.f(idx(idx_max));

    else
        self.para_alpha.kernel.para.mu = self.mufix;
    end

    s_alpha = exp(mean(log(s), 2)) - pred_s_xi(self.f, self.para_xi);
    s_alpha_min = min(s_alpha);
    s_alpha = s_alpha - s_alpha_min;
    opt.maxopt = 'raw';
    % self.para_alpha.kernel.para.h = find_max_around_xt(self.f, s_alpha, self.para_alpha.kernel.para.mu * (1:self.kh)', min(self.width_candidate, self.para_alpha.kernel.para.mu / 2), opt) + s_alpha_min;
    self.para_alpha.kernel.para.h = find_max_around_xt(self.f, s_alpha, get_harmonic_mu_sigma(self.para_alpha.kernel.para.mu, 1, self.ks, self.kh, self.para_alpha.no_zero_peak)', min(self.width_candidate, self.para_alpha.kernel.para.mu / 2), opt) + s_alpha_min;
    self.para_alpha.kernel.para.h = self.para_alpha.kernel.para.h(:);
    self.para_alpha.kernel.para.h(1:self.k0 + self.ks) = 1e-1 * self.para_alpha.kernel.para.h(1:self.k0 + self.ks);
    self.para_alpha.isubharmonic = false;

end

self.para_alpha.kernel.para.h(isnan(self.para_alpha.kernel.para.h) | self.para_alpha.kernel.para.h < 0) = 0;

% fit_alpha(self);

%% remove the effect of trend
% eliminate the xi trend from the height of alpha
% self.para_alpha.kernel.para.h=self.para_alpha.kernel.para.h-pred_s_xi(self.para_alpha.kernel.para.mu*(1:self.kh),self.para_xi).';
% move above
end

%% find first peak
%{
    [self.para_alpha.kernel.para.h(1), self.para_alpha.kernel.para.mu] = findmax_spectrum(log10(self.s - pred_s_xi(self.f, self.para_xi)), self.f, range); %
    %% find other peaks
    % avoid the alpha peak locate larger then the max frequency of
    % sample
    fmax = max(self.f);

    for i = 2:self.kh

        if self.para_alpha.kernel.para.mu * (i) > fmax
            self.kh = i - 1;
            self.para_alpha.kernel.para.h(i:end) = [];
            break;
        end

        range = self.para_alpha.kernel.para.mu * (i) + [-self.width_candidate, self.width_candidate];
        [self.para_alpha.kernel.para.h(i)] = findmax_spectrum(log10(self.s - pred_s_xi(self.f, self.para_xi)), self.f, range);
    end

    % self.para_alpha.kernel.para.h(self.para_alpha.kernel.para.h<0)=0;
    self.para_alpha.kernel.para.h = 10 .^ self.para_alpha.kernel.para.h;

    if k_found<self.k
        self.para_alpha.kernel.para.sigma=self.para_alpha.kernel.para.sigma(1:k_found);
        self.para_alpha.kernel.para.nu=self.para_alpha.kernel.para.nu(1:k_found);
        self.para_alpha.kernel.para.d=self.para_alpha.kernel.para.d(1:k_found);
        self.para_alpha.kernel.para.b=self.para_alpha.kernel.para.b(1:k_found);

        self.para_alpha.kernel.ub.sigma=self.para_alpha.kernel.ub.sigma(1:k_found);
        self.para_alpha.kernel.ub.nu=self.para_alpha.kernel.ub.nu(1:k_found);
        self.para_alpha.kernel.ub.d=self.para_alpha.kernel.ub.d(1:k_found);
        self.para_alpha.kernel.ub.b=self.para_alpha.kernel.ub.b(1:k_found);

        self.para_alpha.kernel.lb.sigma=self.para_alpha.kernel.lb.sigma(1:k_found);
        self.para_alpha.kernel.lb.nu=self.para_alpha.kernel.lb.nu(1:k_found);
        self.para_alpha.kernel.lb.d=self.para_alpha.kernel.lb.d(1:k_found);
        self.para_alpha.kernel.lb.b=self.para_alpha.kernel.lb.b(1:k_found);

    end
%}
