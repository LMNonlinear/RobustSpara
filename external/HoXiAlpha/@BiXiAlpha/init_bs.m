% Author: Ying Wang, Min Li
% Create Time: 2025
% Copyright (c): 2020-2025 Ying Wang, yingwangrigel@gmail.com,
%                Min Li, minli.231314@gmail.com
% Joint China-Cuba LAB, UESTC, Hangzhou Dianzi University
% License: GNU General Public License v3.0 (see LICENSE file)

function init_bs(self)

init_para_bs(self);
bound_bs(self);
bs = mean(self.bs, 3);
bsre = real(bs);
bsim = imag(bs);

% switch self.hos_components
%     case 'xi+alpha'
%         mu1 = [self.para_xi.kernel.para.mu; get_harmonic_mu_sigma(self.para_alpha.kernel.para.mu, 0, self.ks, self.kh, self.para_alpha.no_zero_peak)'];
%         mu2 = get_ndgrid([mu1(:), mu1(:)], 'list');
%     case 'alpha'
%         mu1 = get_harmonic_mu_sigma(self.para_alpha.kernel.para.mu, 0, self.ks, self.kh, self.para_alpha.no_zero_peak)';
%         mu2 = get_ndgrid([mu1(:), mu1(:)], 'list');
%     case 'xi'
%         mu1 = self.para_xi.kernel.para.mu;
%         mu2 = get_ndgrid([mu1(:), mu1(:)], 'list');
%     otherwise
%         error('Unknown HOS components');
% end
if strcmpi(self.peak_relation, 'harmonic')
    mu1 = [self.para_xi.kernel.para.mu; get_harmonic_mu_sigma(self.para_alpha.kernel.para.mu, 0, self.ks, self.kh, self.para_alpha.no_zero_peak)'];
    mu2 = get_ndgrid([mu1(:), mu1(:)], 'list');
else
    mu1=[self.para_xi.kernel.para.mu;self.para_alpha.kernel.para.mu(:)];
    mu2 = get_ndgrid([mu1(:), mu1(:)], 'list');
end
self.para_bs.kernel.para.h(:, :, 1) = reshape(find_max_around_xt(cat(3, self.fx, self.fy), bsre, mu2, self.width_candidate), self.N1, self.N1);
self.para_bs.kernel.para.h(:, :, 2) = reshape(find_max_around_xt(cat(3, self.fx, self.fy), bsim, mu2, self.width_candidate), self.N1, self.N1);

% if strcmpi(self.hos_components, 'xi+alpha')
%     self.para_bs.kernel.para.h(1:self.ks + self.k0, :, :) = 0.1 * self.para_bs.kernel.para.h(1:self.ks + self.k0, :, :);
%     self.para_bs.kernel.para.h(:, 1:self.ks + self.k0, :) = 0.1 * self.para_bs.kernel.para.h(:, 1:self.ks + self.k0, :);
% elseif strcmpi(self.hos_components, 'alpha')
%     self.para_bs.kernel.para.h(1:self.ks + self.k0, :, :) = 0; %* self.para_bs.kernel.para.h(1:self.ks + 1, :, :)
%     self.para_bs.kernel.para.h(:, 1:self.ks + self.k0, :) = 0; %* self.para_bs.kernel.para.h(:, 1:self.ks + 1, :)
% elseif strcmpi(self.hos_components, 'xi')
%     self.para_bs.kernel.para.h(2:end, :, :) = 0; % * self.para_bs.kernel.para.h(self.ks + 2:end, :, :)
%     self.para_bs.kernel.para.h(:, 2:end, :) = 0; % * self.para_bs.kernel.para.h(:, self.ks + 2:end, :)
% else
%     error('Unknown HOS components');
% end
% self.para_bs.kernel.para.h(1:self.ks + self.k0, :, :) = 0.1 * self.para_bs.kernel.para.h(1:self.ks + self.k0, :, :);
% self.para_bs.kernel.para.h(:, 1:self.ks + self.k0, :) = 0.1 * self.para_bs.kernel.para.h(:, 1:self.ks + self.k0, :);

self.para_bs.kernel.para.h(isnan(self.para_bs.kernel.para.h)) = 0;

idx_triu = logical(mat2triu(true(size(self.para_bs.kernel.para.h)), 1, false));
self.para_bs.kernel.para.h(idx_triu) = 0;
end

