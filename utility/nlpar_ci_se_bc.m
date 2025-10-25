function [ci, se] = nlpar_ci_se_bc(bc, resid, jac_r2x, jac_bc2x)
lambda = 1e-3;
alpha = 0.05;

Nr = size(jac_r2x, 1);
Nbc = size(jac_bc2x, 1);
v = Nr - Nbc;
rmse = norm(resid) .^ 2 / v;
% lambda=lambda*trace(jac_r2x'*jac_r2x)/Nr;
% jac_r2x_inv=pinv_ridge(jac_r2x',lambda)';
jac_r2x_inv = pinv_ridge(jac_r2x, lambda);
% lambda=trace(jac_bc2x'*jac_bc2x)/Nbc;
% jac_bc2x_inv=pinv_ridge(jac_bc2x',lambda)';
% J=jac_bc2x*jac_r2x_inv;
% [~, R_bc2x] = qr(jac_bc2x, 0);
% [~, R_r2x_inv] = qr(jac_r2x_inv, 0);
% R_bc2x
% se=jac_bc2x*(jac_r2x_inv*(rmse.*eye(Nr)));
% se=sqrt(rmse*diag(jac_bc2x*(jac_r2x_inv*jac_r2x_inv')*jac_bc2x'));
% tic
se = sqrt(rmse * diag(jac_bc2x * (jac_r2x_inv * jac_r2x_inv') * jac_bc2x'));
% toc

% tic
% se2=sqrt(rmse*sum((jac_bc2x*jac_r2x_inv).^2,2));
% toc

% se=sqrt(diag(jac_bc2x*(jac_r2x_inv*(rmse.*eye(Nr))*jac_r2x_inv')*jac_bc2x'));
% Calculate confidence interval
delta = se * tinv(1 - alpha / 2, v);
ci = [(bc(:) - delta) (bc(:) + delta)];
end

% J'*(sigma*eye(N))*J
