% Author: Ying Wang, Min Li
% Create Time: 2025
% Copyright (c): 2020-2025 Ying Wang, yingwangrigel@gmail.com,
%                Min Li, minli.231314@gmail.com
% Joint China-Cuba LAB, UESTC, Hangzhou Dianzi University
% License: GNU General Public License v3.0 (see LICENSE file)

function CM=interpColor(CList,n)
Ci=1:size(CList,1);
Cq=linspace(1,size(CList,1),n);
CM=[interp1(Ci,CList(:,1),Cq,'pchip')',...
    interp1(Ci,CList(:,2),Cq,'pchip')',...
    interp1(Ci,CList(:,3),Cq,'pchip')'];
end