% Author: Ying Wang, Min Li
% Create Time: 2025
% Copyright (c): 2020-2025 Ying Wang, yingwangrigel@gmail.com,
%                Min Li, minli.231314@gmail.com
% Joint China-Cuba LAB, UESTC, Hangzhou Dianzi University
% License: GNU General Public License v3.0 (see LICENSE file)

function [smax,fmax]=findmax_spectrum(s,f,range)
idx=find(f>=range(1)&f<=range(2));
f=f(idx,:);
s=s(idx,:);
% s_mirror=s;
% s=[s_mirror;s;s_mirror];
% [y(4:-1:2) y y(end-1:-1:end-3)]
% findpeaks(s)
[smax,idx_max]=max(s);
fmax=f(idx_max);


end


