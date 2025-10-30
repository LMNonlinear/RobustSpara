% Author: Ying Wang, Min Li
% Create Time: 2025
% Copyright (c): 2020-2025 Ying Wang, yingwangrigel@gmail.com,
%                Min Li, minli.231314@gmail.com
% Joint China-Cuba LAB, UESTC, Hangzhou Dianzi University
% License: GNU General Public License v3.0 (see LICENSE file)


function [shat,scomp]=pred_s(f,ks,kh,para_xi,para_alpha)
% scomp: different with comp, scomp is with the scale factor, but comp is
%        unit height

if nargout>1
    [s_xi,scomp_xi]=pred_s_xi(f,para_xi);

    [s_alpha,scomp_alpha] = pred_s_alpha(f,ks,kh,para_alpha);%,comp gaussian_kernel(f,mu,sigma).*para_alpha.kernel.alpha.';
    shat=s_xi(:)+s_alpha;
    scomp=[scomp_xi,scomp_alpha];
else
    [s_xi]=pred_s_xi(f,para_xi);
    [s_alpha] = pred_s_alpha(f,ks,kh,para_alpha);%,comp gaussian_kernel(f,mu,sigma).*para_alpha.kernel.alpha.';
    shat=s_xi(:)+s_alpha;
    scomp=[];
end
end




% function [shat,scomp]=pred_s(f,kh,A,B,nu,U,mu,sigma,alpha)
% %         if nargin<2|| isempty(kh)
% %             kh=self.kh;
% %         end
% %         if nargin<3||isempty(mu)   
% %             mu=self.mu;
% %         end
% %         if nargin<4||isempty(sigma)
% %             sigma=self.sigma;
% %         end
% %         if nargin<5||isempty(alpha)
% %             alpha=self.alpha;
% %         end
%         sxi=studentt_waveform(f,A,B,nu,U);

%         [mu,sigma]=get_harmonic_mu_sigma(mu,sigma,kh);
%         shat = sxi(:)+gaussian_kernel(f,mu,sigma)*alpha.';
%         if nargout>1
%             scomp=[sxi,alpha.*gaussian_kernel(f,mu,sigma)];
%         end
%     end