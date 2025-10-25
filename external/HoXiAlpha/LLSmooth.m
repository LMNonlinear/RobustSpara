function [yq, L, dbg] = LLSmooth(x, y, h, xq)
% Local Linear Smoothing (Vectorized for d = 1)
%
% Inputs:
%   x  - Predictor variable, n x 1
%   y  - Response variable, n x 1
%   h  - Bandwidth (scalar or 1 x 1)
%   xq - Query points, nq x 1 (optional, default: x)
%
% Outputs:
%   yq  - Smoothed values at query points, nq x 1
%   L   - Local linear weights, nq x n
%   dbg - Debugging information (sum of weights)

if nargin < 4 || isempty(xq)
    xq = x;
end

% Ensure inputs are column vectors
% x = x(:);
% y = y(:);
% xq = xq(:);

n = length(x);
nq = length(xq);

% Compute scaled differences (nq x n)
D = (xq(:,1) - x') ;  % Broadcasting: nq x n

% Compute kernel weights using Gaussian kernel (nq x n)
W = gaussian_kernel(D,h);   % Already vectorized for d=1

% Compute necessary sums
S0 = sum(W, 2);            % nq x 1
S1 = sum(W .* D, 2);       % nq x 1
S2 = sum(W .* D.^2, 2);    % nq x 1
Sy = sum(W .* y, 2);      % nq x 1
Sxy = sum(W .* D .* y, 2);% nq x 1

% Compute denominator for beta0
denom = S0 .* S2 - S1.^2;  % nq x 1

% Handle potential division by zero
eps_val = 1e-12;
valid = denom > eps_val;

% Initialize yq and L
yq = zeros(nq, 1,size(y,3));
L = zeros(nq, n);

% Compute yq for valid query points
yq(valid,:,:) = (S2(valid) .* Sy(valid,:,:) - S1(valid) .* Sxy(valid,:,:)) ./ denom(valid);

% Fallback to NW smoothing where denominator is too small
fallback = ~valid;
if any(fallback)
    yq(fallback) = Sy(fallback) ./ S0(fallback);
end
if nargout>1
    % Compute Local Linear Weights
    % L = W * (S2 * I - S1 * D') / denom
    % Since d=1, I is scalar 1
    % L = (W .* S2 - W .* D * S1) / denom
    % Expand to nq x n
    L(valid, :) = (W(valid, :) .* S2(valid)) - (W(valid, :) .* D(valid, :) * S1(valid));
    L(valid, :) = L(valid, :) ./ denom(valid);

    % For fallback points, use NW weights
    if any(fallback)
        L(fallback, :) = W(fallback, :) ./ S0(fallback);
    end
end
% Debugging information
dbg.s = S0;

end
