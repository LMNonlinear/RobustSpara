% Author: Ying Wang, Min Li
% Create Time: 2025
% Copyright (c): 2020-2025 Ying Wang, yingwangrigel@gmail.com,
%                Min Li, minli.231314@gmail.com
% Joint China-Cuba LAB, UESTC, Hangzhou Dianzi University
% License: GNU General Public License v3.0 (see LICENSE file)

function [shat,scomp]=pred_s_alpha(f,ks,kh,para_alpha)%,comp
% scomp: different with comp, scomp is with the scale factor, but comp is
%        unit height
% [mu,sigma]=get_harmonic_mu_sigma(para_alpha.kernel.para.mu,para_alpha.kernel.para.sigma,ks,kh);
% comp=gaussian_kernel(f,mu,sigma);

% if nargout>1
%     scomp = comp.*para_alpha.kernel.alpha;
%     shat  = sum(scomp,2);
% else
%     shat = comp*para_alpha.kernel.alpha.';
% end
if nargout>1
    [shat,scomp]=para_alpha.kernel.eval(f);
else
    [shat]=para_alpha.kernel.eval(f);
end
end




