% Author: Ying Wang, Min Li
% Create Time: 2025
% Copyright (c): 2020-2025 Ying Wang, yingwangrigel@gmail.com,
%                Min Li, minli.231314@gmail.com
% Joint China-Cuba LAB, UESTC, Hangzhou Dianzi University
% License: GNU General Public License v3.0 (see LICENSE file)

function y=gaussian_kernel(x,mu,sigma)
if nargin<1
    error(message('gaussian_kernel:TooFewInputs'));
end
if nargin < 2
    mu = 0;
end
if nargin < 3
    sigma = 1;
end

% Return NaN for out of range parameters.
sigma(sigma <= 0) = NaN;

try
    y = exp(-0.5 * ((x - mu)./sigma).^2) ;%./ (sqrt(2*pi) .* sigma)
catch
    error(message('stats:normpdf:InputSizeMismatch'));
end

end