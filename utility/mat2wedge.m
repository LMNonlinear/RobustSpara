% Author: Ying Wang, Min Li
% Create Time: 2025
% Copyright (c): 2020-2025 Ying Wang, yingwangrigel@gmail.com,
%                Min Li, minli.231314@gmail.com
% Joint China-Cuba LAB, UESTC, Hangzhou Dianzi University
% License: GNU General Public License v3.0 (see LICENSE file)

function [matwedge,idx]=mat2wedge(data,isvec,keeptip)
if nargin<3||(keeptip)
    keeptip=true;
end
if nargin<2||isempty(isvec)
    isvec=true;
end
% [matwedge,idx]=mat2wedge(randn(5),false)

[r,c,p]=size(data);

idx=idx_wedge(r,keeptip);

if isvec
    matwedge=data(idx);
else
    matwedge=zeros(r);
    matwedge(idx)=data(idx);
end


end



