% Author: Ying Wang, Min Li
% Create Time: 2025
% Copyright (c): 2020-2025 Ying Wang, yingwangrigel@gmail.com,
%                Min Li, minli.231314@gmail.com
% Joint China-Cuba LAB, UESTC, Hangzhou Dianzi University
% License: GNU General Public License v3.0 (see LICENSE file)

function [lb, ub] = bound_bs(self)

% if nargin<2||isempty(region_bs)
%     region_bs=[];
% end
if strcmpi(self.peak_relation, 'harmonic')
    N1MAX = self.N1MAX;
    N2MAX = self.N2MAX;

    maxbsre = 2.5 * max(abs(real(self.bs)), [], 'all');
    maxbsim = 2.5 * max(abs(imag(self.bs)), [], 'all');

    % clb = 1; %;1;0.75
    % cub = 1; %1; 1.25
    clb = 0.75; %;1;0.75
    cub = 1.25; %1; 1.25

    self.para_bs.kernel.lb.h = cat(3, -maxbsre * ones(N1MAX, N1MAX), -maxbsim * ones(N1MAX, N1MAX));
    self.para_bs.kernel.lb.mu = [self.para_xi.kernel.para.mu; self.para_alpha.kernel.para.mu];
    self.para_bs.kernel.lb.sigma = [self.para_xi.kernel.para.sigma, clb; ...
                                        self.para_alpha.kernel.para.sigma, clb];
    self.para_bs.kernel.lb.nu = [self.para_xi.kernel.para.nu; ...
                                     self.para_alpha.kernel.para.nu];
    self.para_bs.kernel.lb.d = [self.para_xi.kernel.para.d; ...
                                    self.para_alpha.kernel.para.d];
    self.para_bs.kernel.lb.b = [self.para_xi.kernel.para.b; ...
                                    self.para_alpha.kernel.para.b];

    self.para_bs.kernel.ub.h = cat(3, maxbsre * ones(N1MAX, N1MAX), maxbsim * ones(N1MAX, N1MAX));
    self.para_bs.kernel.ub.mu = [self.para_xi.kernel.para.mu; self.para_alpha.kernel.para.mu];
    self.para_bs.kernel.ub.sigma = [self.para_xi.kernel.para.sigma, cub; ...
                                        self.para_alpha.kernel.para.sigma, cub];
    self.para_bs.kernel.ub.nu = [self.para_xi.kernel.para.nu; ...
                                     self.para_alpha.kernel.para.nu];
    self.para_bs.kernel.ub.d = [self.para_xi.kernel.para.d; ...
                                    self.para_alpha.kernel.para.d];
    self.para_bs.kernel.ub.b = [self.para_xi.kernel.para.b; ...
                                    self.para_alpha.kernel.para.b];

    if self.para_bs.use_bssigma

        self.para_bs.kernel.lb.sigma = [self.para_xi.kernel.para.sigma, 1; ...
                                            self.para_alpha.kernel.para.sigma ./ 2, 0.75; ];
        self.para_bs.kernel.ub.sigma = [self.para_xi.kernel.para.sigma, 1; ...
                                            self.para_alpha.kernel.para.sigma * 2, 1.25]; %+ 1
    end

    if self.para_bs.use_bsnu

        self.para_bs.kernel.lb.nu = [self.para_xi.kernel.para.nu; ...
                                         self.para_alpha.kernel.lb.nu];
        self.para_bs.kernel.ub.nu = [self.para_xi.kernel.para.nu; ...
                                         self.para_alpha.kernel.ub.nu];

    end

    if self.para_bs.use_bsd

        self.para_bs.kernel.lb.d = [self.para_xi.kernel.para.d; ...
                                        self.para_alpha.kernel.lb.d];
        self.para_bs.kernel.ub.d = [self.para_xi.kernel.para.d; ...
                                        self.para_alpha.kernel.ub.d];

    end

    if self.para_bs.use_bsb

        self.para_bs.kernel.lb.b = [self.para_xi.kernel.lb.b; ...
                                        self.para_alpha.kernel.lb.b];
        self.para_bs.kernel.ub.b = [self.para_xi.kernel.ub.b; ...
                                        self.para_alpha.kernel.ub.b];

    end

    switch self.hos_components
        case 'xi+alpha'

        case 'alpha'
            self.para_bs.kernel.lb.h(1, :, :) = 0;
            self.para_bs.kernel.lb.h(:, 1, :) = 0;
            self.para_bs.kernel.ub.h(1, :, :) = 0;
            self.para_bs.kernel.ub.h(:, 1, :) = 0;
        case 'xi'
            self.para_bs.kernel.lb.h(2:end, :, :) = 0;
            self.para_bs.kernel.lb.h(:, 2:end, :) = 0;
            self.para_bs.kernel.ub.h(2:end, :, :) = 0;
            self.para_bs.kernel.ub.h(:, 2:end, :) = 0;
        otherwise
            error('no such method')
    end

    % reduce the number of parameters by the property of the bispectrum, only
    % keep lower half
    idx_triu = logical(mat2triu(true(N1MAX, N1MAX, 2), 1, false));
    self.para_bs.kernel.lb.h(idx_triu) = 0;
    self.para_bs.kernel.ub.h(idx_triu) = 0;

    if self.para_bs.no_xi_alpha
        idx_xialpha = mat2triu(false(N1MAX, N1MAX, 2), 1, false);
        idx_xialpha(1, 2:end, :) = true;
        idx_xialpha(2:end, 1, :) = true;
        idx_xialpha = logical(idx_xialpha);
        self.para_bs.kernel.lb.h(idx_xialpha) = 0;
        self.para_bs.kernel.ub.h(idx_xialpha) = 0;
    end

    %
    if self.para_bs.no_subharmonic_off_diag
        idx_subharmonic_offdiag = false(N1MAX, N1MAX, 2); %mat2triu(,1);
        idx_subharmonic_offdiag(2:self.k0 + self.ksmax + 1, :, :) = true;
        idx_subharmonic_offdiag(:, 2:self.k0 + self.ksmax + 1, :) = true;
        idx_subharmonic_offdiag(logical(diag2full(ones(size(idx_subharmonic_offdiag, 1), 2)))) = false;
        self.para_bs.kernel.lb.h(idx_subharmonic_offdiag) = 0;
        self.para_bs.kernel.ub.h(idx_subharmonic_offdiag) = 0;
    end

    %
    mus = get_harmonic_mu_sigma(self.para_alpha.kernel.para.mu, 0, self.ks, self.kh, self.para_alpha.no_zero_peak); %, self.df

    switch self.hos_components
        case 'xi+alpha'
            mus = [0, mus];
        case 'alpha'
            % mus = mus;
        case 'xi'

        otherwise
            error('no such method')
    end

    [mux, muy] = meshgrid(mus, mus);
    mask = logical(sum(abs(self.para_bs.kernel.lb.h) ~= 0, 3));

    mux(~mask) = 0;
    muy(~mask) = 0;
    mu = cat(3, mux, muy);
    mu_id = round(mu ./ self.para_alpha.kernel.para.mu);
    idx_outside = sum(mu_id, 3) > self.kh;
    idx_outside = cat(3, idx_outside, idx_outside);
    self.para_bs.kernel.lb.h(idx_outside) = 0;
    self.para_bs.kernel.ub.h(idx_outside) = 0;

    % assign zero to the peak at zero
    if ~self.para_alpha.no_zero_peak
        self.para_bs.kernel.lb.h(2, 2, :) = 0;
        self.para_bs.kernel.ub.h(2, 2, :) = 0;
    end

    % assign the fixed parameters
    idx_fix_pre = self.para_bs.kernel.lb.h == self.para_bs.kernel.ub.h;
    self.para_bs.kernel.para.h(idx_fix_pre) = self.para_bs.kernel.lb.h(idx_fix_pre);

    switch self.region_bs
        case {'', 'all'}
            idx_fix = false(N1MAX, N1MAX, 2);
        case 'mu+diag'
            idx_fix = logical(~(eye(N1MAX) | [zeros(N1MAX, self.k0 + self.ks), ones(N1MAX, 1), zeros(N1MAX, N1MAX - 1 - self.k0 - self.ks)]) .* ones(1, 1, 2));
            idx_fix = idx_fix & permute(idx_fix, [2, 1, 3]);
        case 'boundary'
            idx_fix = logical([zeros(N1MAX, 1), ones(N1MAX, N1MAX - 1)] .* ones(1, 1, 2));
            idx_fix = idx_fix & permute(idx_fix, [2, 1, 3]);
        case 'diag'
            idx_fix = (~eye(N1MAX) .* ones(1, 1, 2));
        case 'offdiag'
            idx_fix = (eye(N1MAX) .* ones(1, 1, 2));
        otherwise
            error('no such region_bs')
    end

    idx_fix = logical(idx_fix);
    % idx_fix = idx_fix| idx_fix_pre;%
    self.para_bs.kernel.lb.h(idx_fix) = self.para_bs.kernel.para.h(idx_fix);
    self.para_bs.kernel.ub.h(idx_fix) = self.para_bs.kernel.para.h(idx_fix);

    if ~isempty(self.para_bs.kernel.fix)
        self.para_bs.kernel.vlb(~isnan(self.para_bs.kernel.vfix)) = self.para_bs.kernel.vfix(~isnan(self.para_bs.kernel.vfix));
        self.para_bs.kernel.vub(~isnan(self.para_bs.kernel.vfix)) = self.para_bs.kernel.vfix(~isnan(self.para_bs.kernel.vfix));
    end

    lb = self.para_bs.kernel.vlb;
    ub = self.para_bs.kernel.vub;

elseif strcmpi(self.peak_relation, 'free')
    N1MAX = self.N1MAX;
    N2MAX = self.N2MAX;

    maxbsre = 2.5 * max(abs(real(self.bs)), [], 'all');
    maxbsim = 2.5 * max(abs(imag(self.bs)), [], 'all');

    clb = 1; %;1;0.75
    cub = 1; %1; 1.25
    % clb = 0.75; %;1;0.75
    % cub = 1.25; %1; 1.25

    self.para_bs.kernel.lb.h = cat(3, -maxbsre * ones(N1MAX, N1MAX), -maxbsim * ones(N1MAX, N1MAX));
    self.para_bs.kernel.lb.mu = [self.para_xi.kernel.para.mu; self.para_alpha.kernel.para.mu];

    % self.para_bs.kernel.lb.sigma = [self.para_xi.kernel.para.sigma, -c; ...
    %                                   self.para_alpha.kernel.para.sigma, -c];
    % self.para_bs.kernel.lb.sigma = [self.para_xi.kernel.para.sigma, clb; ...
    %     self.para_alpha.kernel.para.sigma, clb];
    % sigma_lb1=[self.para_xi.kernel.para.sigma;self.para_alpha.kernel.para.sigma];
    % self.para_bs.kernel.lb.sigma = get_ndgrid([sigma_lb1,sigma_lb1], 'list');

    % self.para_bs.kernel.ub.h = cat(3, maxbsre * ones(N1MAX, N1MAX), maxbsim * ones(N1MAX, N1MAX));
    % self.para_bs.kernel.ub.mu = [self.para_xi.kernel.para.mu; self.para_alpha.kernel.para.mu];
    mu1 = [self.para_xi.kernel.para.mu; self.para_alpha.kernel.para.mu];
    self.para_bs.kernel.lb.mu = get_ndgrid([mu1, mu1], 'list');

    sigma_lb1 = [self.para_xi.kernel.para.sigma; self.para_alpha.kernel.para.sigma];
    self.para_bs.kernel.lb.sigma = get_ndgrid([sigma_lb1, sigma_lb1], 'list');

    nu_lb1 = [self.para_xi.kernel.para.nu; self.para_alpha.kernel.para.nu];
    self.para_bs.kernel.lb.nu = get_ndgrid([nu_lb1, nu_lb1], 'list');

    d_lb1 = [self.para_xi.kernel.para.d; self.para_alpha.kernel.para.d];
    self.para_bs.kernel.lb.d = get_ndgrid([d_lb1, d_lb1], 'list');

    b_lb1 = [self.para_xi.kernel.para.b; self.para_alpha.kernel.para.b];
    self.para_bs.kernel.lb.b = get_ndgrid([b_lb1, b_lb1], 'list');

    self.para_bs.kernel.ub.h = cat(3, maxbsre * ones(N1MAX, N1MAX), maxbsim * ones(N1MAX, N1MAX));
    % self.para_bs.kernel.ub.mu = [self.para_xi.kernel.para.mu; self.para_alpha.kernel.para.mu];
    mu1 = [self.para_xi.kernel.para.mu; self.para_alpha.kernel.para.mu];
    self.para_bs.kernel.ub.mu = get_ndgrid([mu1, mu1], 'list');

    sigma_ub1 = [self.para_xi.kernel.para.sigma; self.para_alpha.kernel.para.sigma];
    self.para_bs.kernel.ub.sigma = get_ndgrid([sigma_ub1, sigma_ub1], 'list');

    nu_ub1 = [self.para_xi.kernel.para.nu; self.para_alpha.kernel.para.nu];
    self.para_bs.kernel.ub.nu = get_ndgrid([nu_ub1, nu_ub1], 'list');

    d_ub1 = [self.para_xi.kernel.para.d; self.para_alpha.kernel.para.d];
    self.para_bs.kernel.ub.d = get_ndgrid([d_ub1, d_ub1], 'list');

    b_ub1 = [self.para_xi.kernel.para.b; self.para_alpha.kernel.para.b];
    self.para_bs.kernel.ub.b = get_ndgrid([b_ub1, b_ub1], 'list');

    if self.para_bs.use_bssigma
        %{
    clb = 0.75; %;0.5;
    cub = 1.25; %2.5;
    sigma_lb_xi = max(1/2 * self.para_xi.kernel.para.sigma, 50); %self.para_xi.kernel.lb.sigma;
    sigma_lb_alpha = max(1/2 * self.para_alpha.kernel.para.sigma, 0.5); %self.para_alpha.kernel.lb.sigma; %
    % self.para_bs.kernel.lb.sigma = [sigma_lb_xi, -c; ...
    %                                   sigma_lb_alpha, -c];
    self.para_bs.kernel.lb.sigma = [sigma_lb_xi, 1; ...
                                        sigma_lb_alpha, 1];
    sigma_ub_xi = max(max(max(self.f), 100), 2 * self.para_xi.kernel.para.sigma); %self.para_xi.kernel.ub.sigma; %
    sigma_ub_alpha = max(2 * self.para_alpha.kernel.para.sigma, 2); %self.para_alpha.kernel.ub.sigma; %
    self.para_bs.kernel.ub.sigma = [sigma_ub_xi, 1; ...
                                        sigma_ub_alpha, 1];
    % self.para_bs.kernel.ub.sigma=[sigma_ub_xi,c;...
    %     sigma_ub_alpha,c];
        %}

        % sigma_lb_xi = max(1/2 * self.para_xi.kernel.para.sigma, 50); %self.para_xi.kernel.lb.sigma;
        % sigma_lb_alpha = max(1/2 * self.para_alpha.kernel.para.sigma, 0.5); %self.para_alpha.kernel.lb.sigma; %
        % sigma_ub_xi = max(max(max(self.f), 100), 2 * self.para_xi.kernel.para.sigma); %self.para_xi.kernel.ub.sigma; %
        % sigma_ub_alpha = max(2 * self.para_alpha.kernel.para.sigma, 2); %self.para_alpha.kernel.ub.sigma; %
        % self.para_bs.kernel.lb.sigma = [sigma_lb_xi, 1; ...
        %                                     sigma_lb_alpha, 1];
        % self.para_bs.kernel.ub.sigma = [sigma_ub_xi, 1; ...
        %                                     sigma_ub_alpha, 1];
        % self.para_bs.kernel.lb.sigma = [self.para_xi.kernel.para.sigma, 1; ...
        %     self.para_alpha.kernel.para.sigma ./ 2, 0.75; ];
        % self.para_bs.kernel.ub.sigma = [self.para_xi.kernel.para.sigma, 1; ...
        %     self.para_alpha.kernel.para.sigma, 1.25]; %+ 1

        sigma_lb1 = [self.para_xi.kernel.lb.sigma; self.para_alpha.kernel.lb.sigma];
        self.para_bs.kernel.lb.sigma = get_ndgrid([sigma_lb1, sigma_lb1], 'list');
        sigma_ub1 = [self.para_xi.kernel.ub.sigma; self.para_alpha.kernel.ub.sigma];
        self.para_bs.kernel.ub.sigma = get_ndgrid([sigma_ub1, sigma_ub1], 'list');
    end

    if self.para_bs.use_bsnu
        %{
    clb = 0.75; %;0.5;
    cub = 1.25; %2.5;
    self.para_bs.kernel.lb.nu = [clb * self.para_xi.kernel.para.nu; ...
                                     clb * self.para_alpha.kernel.para.nu];
    self.para_bs.kernel.ub.nu = [cub * self.para_xi.kernel.para.nu; ...
                                     cub * self.para_alpha.kernel.para.nu];
        %}
        % self.para_bs.kernel.lb.nu = [self.para_xi.kernel.para.nu; ...
        %     self.para_alpha.kernel.lb.nu];
        % self.para_bs.kernel.ub.nu = [self.para_xi.kernel.para.nu; ...
        %     self.para_alpha.kernel.ub.nu];
        % self.para_bs.kernel.lb.nu = [self.para_xi.kernel.para.nu; ...
        %                                  self.para_alpha.kernel.para.nu];
        % self.para_bs.kernel.ub.nu = [self.para_xi.kernel.para.nu; ...
        %                                  self.para_alpha.kernel.para.nu * 2]; %+ 1
        nu_lb1 = [self.para_xi.kernel.lb.nu; self.para_alpha.kernel.lb.nu];
        self.para_bs.kernel.lb.nu = get_ndgrid([nu_lb1, nu_lb1], 'list');
        nu_ub1 = [self.para_xi.kernel.ub.nu; self.para_alpha.kernel.ub.nu];
        self.para_bs.kernel.ub.nu = get_ndgrid([nu_ub1, nu_ub1], 'list');

    end

    if self.para_bs.use_bsd
        %{
    clb = 0.75; %;0.5;
    cub = 1.25; %2.5;
    self.para_bs.kernel.lb.d = [clb * self.para_xi.kernel.para.d; ...
                                   clb * self.para_alpha.kernel.para.d];
    self.para_bs.kernel.ub.d = [cub * self.para_xi.kernel.para.d; ...
                                   cub * self.para_alpha.kernel.para.d];
        %}
        % self.para_bs.kernel.lb.d = [self.para_xi.kernel.para.d; ...
        %     self.para_alpha.kernel.lb.d];
        % self.para_bs.kernel.ub.d = [self.para_xi.kernel.para.d; ...
        %     self.para_alpha.kernel.ub.d];
        d_lb1 = [self.para_xi.kernel.lb.d; self.para_alpha.kernel.lb.d];
        self.para_bs.kernel.lb.d = get_ndgrid([d_lb1, d_lb1], 'list');
        d_ub1 = [self.para_xi.kernel.ub.d; self.para_alpha.kernel.ub.d];
        self.para_bs.kernel.ub.d = get_ndgrid([d_ub1, d_ub1], 'list');

    end

    if self.para_bs.use_bsb
        %{
    clb = 0.75; %;0.5;
    cub = 1.25; %2.5;
    self.para_bs.kernel.lb.b = [clb * self.para_xi.kernel.para.b; ...
                                   clb * self.para_alpha.kernel.para.b];
    self.para_bs.kernel.ub.b = [cub * self.para_xi.kernel.para.b; ...
                                   cub * self.para_alpha.kernel.para.b];
        %}
        % self.para_bs.kernel.lb.b = [self.para_xi.kernel.lb.b; ...
        %     self.para_alpha.kernel.lb.b];
        % self.para_bs.kernel.ub.b = [self.para_xi.kernel.ub.b; ...
        %     self.para_alpha.kernel.ub.b];
        b_lb1 = [self.para_xi.kernel.lb.b; self.para_alpha.kernel.lb.b];
        self.para_bs.kernel.lb.b = get_ndgrid([b_lb1, b_lb1], 'list');

        b_ub1 = [self.para_xi.kernel.ub.b; self.para_alpha.kernel.ub.b];
        self.para_bs.kernel.ub.b = get_ndgrid([b_ub1, b_ub1], 'list');

    end

    switch self.hos_components
        case 'xi+alpha'

        case 'alpha'
            self.para_bs.kernel.lb.h(1, :, :) = 0;
            self.para_bs.kernel.lb.h(:, 1, :) = 0;
            self.para_bs.kernel.ub.h(1, :, :) = 0;
            self.para_bs.kernel.ub.h(:, 1, :) = 0;
        case 'xi'
            self.para_bs.kernel.lb.h(2:end, :, :) = 0;
            self.para_bs.kernel.lb.h(:, 2:end, :) = 0;
            self.para_bs.kernel.ub.h(2:end, :, :) = 0;
            self.para_bs.kernel.ub.h(:, 2:end, :) = 0;
        otherwise
            error('no such method')
    end

    % reduce the number of parameters by the property of the bispectrum, only
    % keep lower half
    idx_triu = logical(mat2triu(true(N1MAX, N1MAX, 2), 1, false));
    self.para_bs.kernel.lb.h(idx_triu) = 0;
    self.para_bs.kernel.ub.h(idx_triu) = 0;

    if self.para_bs.no_xi_alpha
        idx_xialpha = mat2triu(false(N1MAX, N1MAX, 2), 1, false);
        idx_xialpha(1, 2:end, :) = true;
        idx_xialpha(2:end, 1, :) = true;
        idx_xialpha = logical(idx_xialpha);
        self.para_bs.kernel.lb.h(idx_xialpha) = 0;
        self.para_bs.kernel.ub.h(idx_xialpha) = 0;
    end

    %
    if self.para_bs.no_subharmonic_off_diag
        idx_subharmonic_offdiag = false(N1MAX, N1MAX, 2); %mat2triu(,1);
        idx_subharmonic_offdiag(2:self.k0 + self.ksmax + 1, :, :) = true;
        idx_subharmonic_offdiag(:, 2:self.k0 + self.ksmax + 1, :) = true;
        idx_subharmonic_offdiag(logical(diag2full(ones(size(idx_subharmonic_offdiag, 1), 2)))) = false;
        self.para_bs.kernel.lb.h(idx_subharmonic_offdiag) = 0;
        self.para_bs.kernel.ub.h(idx_subharmonic_offdiag) = 0;
    end

    %{
    mus = get_harmonic_mu_sigma(self.para_alpha.kernel.para.mu, 0, self.ks, self.kh, self.para_alpha.no_zero_peak); %, self.df

    switch self.hos_components
        case 'xi+alpha'
            mus = [0, mus];
        case 'alpha'
            % mus = mus;
        case 'xi'

        otherwise
            error('no such method')
    end

    [mux, muy] = meshgrid(mus, mus);
    mask = logical(sum(abs(self.para_bs.kernel.lb.h) ~= 0, 3));

    mux(~mask) = 0;
    muy(~mask) = 0;
    mu = cat(3, mux, muy);
    mu_id = round(mu ./ self.para_alpha.kernel.para.mu);
    idx_outside = sum(mu_id, 3) > self.kh;
    idx_outside = cat(3, idx_outside, idx_outside);
    self.para_bs.kernel.lb.h(idx_outside) = 0;
    self.para_bs.kernel.ub.h(idx_outside) = 0;
    %}
    % assign zero to the peak at zero
    if ~self.para_alpha.no_zero_peak
        self.para_bs.kernel.lb.h(2, 2, :) = 0;
        self.para_bs.kernel.ub.h(2, 2, :) = 0;
    end

    % assign the fixed parameters
    idx_fix_pre = self.para_bs.kernel.lb.h == self.para_bs.kernel.ub.h;
    self.para_bs.kernel.para.h(idx_fix_pre) = self.para_bs.kernel.lb.h(idx_fix_pre);

    switch self.region_bs
        case {'', 'all'}
            idx_fix = false(N1MAX, N1MAX, 2);
        case 'mu+diag'
            idx_fix = logical(~(eye(N1MAX) | [zeros(N1MAX, self.k0 + self.ks), ones(N1MAX, 1), zeros(N1MAX, N1MAX - 1 - self.k0 - self.ks)]) .* ones(1, 1, 2));
            idx_fix = idx_fix & permute(idx_fix, [2, 1, 3]);
        case 'boundary'
            idx_fix = logical([zeros(N1MAX, 1), ones(N1MAX, N1MAX - 1)] .* ones(1, 1, 2));
            idx_fix = idx_fix & permute(idx_fix, [2, 1, 3]);
        case 'diag'
            idx_fix = (~eye(N1MAX) .* ones(1, 1, 2));
        case 'offdiag'
            idx_fix = (eye(N1MAX) .* ones(1, 1, 2));
        otherwise
            error('no such region_bs')
    end

    idx_fix = logical(idx_fix);
    % idx_fix = idx_fix| idx_fix_pre;%
    self.para_bs.kernel.lb.h(idx_fix) = self.para_bs.kernel.para.h(idx_fix);
    self.para_bs.kernel.ub.h(idx_fix) = self.para_bs.kernel.para.h(idx_fix);

    if ~isempty(self.para_bs.kernel.fix)
        self.para_bs.kernel.vlb(~isnan(self.para_bs.kernel.vfix)) = self.para_bs.kernel.vfix(~isnan(self.para_bs.kernel.vfix));
        self.para_bs.kernel.vub(~isnan(self.para_bs.kernel.vfix)) = self.para_bs.kernel.vfix(~isnan(self.para_bs.kernel.vfix));
    end

    lb = self.para_bs.kernel.vlb;
    ub = self.para_bs.kernel.vub;

end

end

