% Author: Ying Wang, Min Li
% Create Time: 2025
% Copyright (c): 2020-2025 Ying Wang, yingwangrigel@gmail.com,
%                Min Li, minli.231314@gmail.com
% Joint China-Cuba LAB, UESTC, Hangzhou Dianzi University
% License: GNU General Public License v3.0 (see LICENSE file)

function [y, ycomp] = kernel_harmonic_t_flattop(x, h, mu, sigma, nu, d, b)
% D = bs_dist(x, mu(imu, :), Sigma(:, :, idsigma), nu(idsigma));
% p = size(x, 2);
% sigma = sigma(:, 1, :);
% Sigma=diag([sigma(1);sigma(1);sigma(2)]);
p = 2;

if size(sigma, 2) == 3
    % sigma = permute([sigma(1, :, :); sigma(1, :, :); sigma(2, :, :)], [2, 1, 3]); %
    x = [x, sum(x, 2)]; %x,y,x+y
    mu = [mu, sum(mu, 2)];
    % p = 3;
    p = 1;
% elseif size(sigma, 1) > 2
%     error('unsupport dim > 2')
end

sigma = permute(sigma, [3, 1, 2]); % Nmu*Nd*1->1*Nmu*Nd
mu = permute(mu, [3, 1, 2]); % Nmu*Nd*1->1*Nmu*Nd
nu = permute(nu, [3, 1, 2]); % Nmu*Nd*1->1*Nmu*Nd
d = permute(d, [3, 1, 2]); % Nmu*Nd*1->1*Nmu*Nd
b = permute(b, [3, 1, 2]); % Nmu*Nd*1->1*Nmu*Nd

x = permute(x, [1, 3, 2]); % Nx*Nd->Nx*1*Nd

% sigma = scale_sigma(sigma, h, nu, d, p);
% h = h .^ (1 ./ p);

if nargout > 1
    % ycomp = h .* (prod((1 + abs((x - mu) ./ sigma) .^ 2) .^ (-1 ./ nu), 3));
    % ycomp = h .* (prod((1 + abs((x - mu) ./ sigma) .^ d) .^ (-1 ./ nu), 3));
    % ycomp = h .* (prod((1 + abs((x - mu) ./ sigma) .^ d) .^ (-1 ./ (nu)), 3));
    % ycomp = h .* (prod((1 + abs((x - mu) ./ sigma) .^ d) .^ (-1 ./ nu), 3) .^ p);
    % ycomp = h .* (prod((1 + abs((x - mu) ./ sigma) .^ 2) .^ (-1 ./ nu), 3) .^ p);
    ycomp = h .* (prod((1 + abs((x - mu) ./ sigma) .^ d) .^ (- nu), 3) .^ p);
    % ycomp = h .* prod((10 .^ -b) ./ (exp((abs((x - mu) ./ sigma)) .^ (1 ./ (nu)))) + 1 ./ (1 + abs((x - mu) ./ sigma) .^ (d)), 3);
    % ycomp = h .* prod((10 .^ -b) ./ (exp(((nu)) .^ (abs((x - mu) ./ sigma)))) + 1 ./ (1 + abs((x - mu) ./ sigma) .^ (d)), 3);
    y = sum(ycomp, 2);

    % if p == 2
    %     ycomp = h .* tkernel(x, mu, sigma, nu, d);
    %     y = sum(ycomp, 2) .^ 2;
    %     ycomp = ycomp .^ 2;
    % elseif p == 3
    %     ycomp = h .* tkernel(x, mu, sigma, nu, d);
    %     y = prod(sum(ycomp, 2), 3);
    % end

else
    % y = sum(h .* prod((1 + abs((x - mu) ./ sigma) .^ 2) .^ (-1 ./ nu), 3), 2);
    % y = sum(h .* prod((1 + abs((x - mu) ./ sigma) .^ d) .^ (-1 ./ nu), 3), 2);
    % y = sum(h .* prod((1 + abs((x - mu) ./ sigma) .^ d) .^ (-1 ./ (nu)), 3), 2);
    % y = sum(h .* prod((1 + abs((x - mu) ./ sigma) .^ d) .^ (-1 ./ nu), 3) .^ p, 2);
    % y = sum(h .* prod((1 + abs((x - mu) ./ sigma) .^ 2) .^ (-1 ./ nu), 3) .^ p, 2);
    y = sum(h .* prod((1 + abs((x - mu) ./ sigma) .^ d) .^ (- nu), 3) .^ p, 2);
    % y = sum(h .* prod((10 .^ -b) ./ (exp((abs((x - mu) ./ sigma)) .^ (1 ./ (nu)))) + 1 ./ (1 + abs((x - mu) ./ sigma) .^ (d)), 3), 2);

    % if p == 2
    %     y = sum(h .* tkernel(x, mu, sigma, nu, d), 2) .^ 2;
    % elseif p == 3
    %     h = permute(h, [1, 2, 4, 3]);
    %     y = prod(sum(h .* tkernel(x, mu, sigma, nu, d), 2), 3);
    %     y = permute(y, [1, 2, 4, 3]);
    % end

end

% if ndims(y)>3
%     y = permute(y, [1, 2, 4, 3]);
% end
end

function sigma = scale_sigma(sigma, h, nu, d, p)
% decouple sigma and nu
sigma = sigma ./ ((2) .^ (1 ./ (nu)) - 1) .^ (1 ./ d);

end

% function sigma = scale_sigma(sigma,h, nu, d, p)
% if p==1
%     h = permute(h, [1, 2, 4, 3]);
%     h = h .^ (1 ./ 3);
% elseif p==2
%     h = h .^ (1 ./ 2);
% else
%     error('unkown p')
% end
%
%
% % decouple sigma and nu
% sigma = sigma ./ ( h.^ (1 ./ (2*nu)) - 1 ) .^ (1 ./ d);
%
% end

function y = tkernel(x, mu, sigma, nu, d)
y = (1 + abs((x - mu) ./ sigma) .^ d) .^ (-1 ./ nu);
end
