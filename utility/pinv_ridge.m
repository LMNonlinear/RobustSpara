% Author: Ying Wang, Min Li
% Create Time: 2025
% Copyright (c): 2020-2025 Ying Wang, yingwangrigel@gmail.com,
%                Min Li, minli.231314@gmail.com
% Joint China-Cuba LAB, UESTC, Hangzhou Dianzi University
% License: GNU General Public License v3.0 (see LICENSE file)

function Jinv=pinv_ridge(J,lambda)
% pinv_ridge=@(J,lambda) (J'*J+lambda*eye(size(J,2)))\J';
% Jinv=(J'*J+lambda*eye(size(J,2)))\J';
% 
N=size(J,1);
% lambda=lambda*trace(J'*J)/N;
lambda=lambda*sum(J.*J,'all')/N;
Jinv=(J'*J+lambda*eye(size(J,2)))\J';
end