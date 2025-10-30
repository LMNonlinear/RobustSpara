% Author: Ying Wang, Min Li
% Create Time: 2025
% Copyright (c): 2020-2025 Ying Wang, yingwangrigel@gmail.com,
%                Min Li, minli.231314@gmail.com
% Joint China-Cuba LAB, UESTC, Hangzhou Dianzi University
% License: GNU General Public License v3.0 (see LICENSE file)

function surfbs(fx, fy, bs,issymcolor)
if nargin<4 || isempty(issymcolor)
    issymcolor=true;
end
p = surf(fx, fy, bs);
view(0, 90)
xlim([min(fx, [], 'all'), max(fx, [], 'all')])
ylim([min(fy, [], 'all'), max(fy, [], 'all')])
set(p, 'edgecolor', 'none')
% shading interp
% grid off
ax=gca;
colormap(ax,turbo) % colormap(flip(hot)) %hot
% colormap(flip(hot))

colorbar
if issymcolor
    maxbs=max(abs(bs),[],"all");
%     clim([-maxbs,maxbs])
end
title('bispectrum')
box off
end
