% Author: Ying Wang, Min Li
% Create Time: 2025
% Copyright (c): 2020-2025 Ying Wang, yingwangrigel@gmail.com,
%                Min Li, minli.231314@gmail.com
% Joint China-Cuba LAB, UESTC, Hangzhou Dianzi University
% License: GNU General Public License v3.0 (see LICENSE file)

function [mu, sigma, order] = get_harmonic_mu_sigma(mu, sigma, ks, kh, no_zero_peak)
% warning('odd subharmonic is not visible in the bispectrum')

order = [1 ./ ((ks:-1:1) + 1), (1:kh)];

mu = mu(1) .* order;

% if ~no_zero_peak
%     mu = [df, mu];
% end
if ~no_zero_peak
    mu = [0, mu];
end

if no_zero_peak
    sigma = [sigma(1) * ones(1, ks + kh)];
else
    sigma = [sigma(1) * ones(1, 1 + ks + kh)];
end

end
