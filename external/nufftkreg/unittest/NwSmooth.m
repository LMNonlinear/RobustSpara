% Author: Ying Wang, Min Li
% Create Time: 2025
% Copyright (c): 2020-2025 Ying Wang, yingwangrigel@gmail.com,
%                Min Li, minli.231314@gmail.com
% Joint China-Cuba LAB, UESTC, Hangzhou Dianzi University
% License: GNU General Public License v3.0 (see LICENSE file)

function [yq,L,dbg]=NwSmooth(x,y,h,xq)
[x,x_mean,x_std]=zscore(x);
if  nargin<4||isempty(xq)
    xq=x;
end
xq=(xq-x_mean)./x_std;
% x=single(x);
% xq=single(xq);
% y=single(y);

[n,dx]=size(x);
[nq,~]=size(xq);
if ~isvector(xq)
    xq=reshape(xq,[nq,1,dx]);
    x=reshape(x,[n,1,dx]);
end
h=permute(h,[dx+2:-1:1]);

D=x-permute(xq,[2,1,3]);
% Dkn=squeeze(gauss_kernel(D,h));
% Dkn=squeeze(epan_kernel(D,h));
Dkn=squeeze(gaussian_kernel(D,h));
% Dx=ones(n,1);
% e1=[1;zeros(n-1,1)];
% L=e1*(Dx'*Dkn*Dx./Dx'*Dkn);
% yq=L*y;
yq=squeeze(sum(Dkn.*y,1)./sum(Dkn,1));
yq=yq';
dbg.s=sum(Dkn,1);
if nargout>1
    L=Dkn./sum(Dkn,1).';%
%       L=diag(Dkn)./sum(Dkn,1).';
end
end


function u=epan_kernel(u,b)
u=u./b;
% u=max(eps,3/4*(1-sum(u.^2,3)));
u=max(0,3/4*(1-sum(u.^2,3)));
end

% function w = gauss_kernel(u,b)
% u=u./b;
% w = exp(-0.5*sum(u.^2,3));
% end


function u=gaussian_kernel(u,b)
u=u./b;
u=(1/sqrt(2*pi))*exp(-0.5*sum((u.^2),3));
end


