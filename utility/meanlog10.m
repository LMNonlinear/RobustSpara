% Author: Ying Wang, Min Li
% Create Time: 2025
% Copyright (c): 2020-2025 Ying Wang, yingwangrigel@gmail.com,
%                Min Li, minli.231314@gmail.com
% Joint China-Cuba LAB, UESTC, Hangzhou Dianzi University
% License: GNU General Public License v3.0 (see LICENSE file)

function mu=meanlog10(x,dim,varargin)
if nargin<2 || isempty(dim)
    dim=find(size(x)>1);
end

mu=10.^mean(log10(x),dim,varargin{:});


end

