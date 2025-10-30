% Author: Ying Wang, Min Li
% Create Time: 2025
% Copyright (c): 2020-2025 Ying Wang, yingwangrigel@gmail.com,
%                Min Li, minli.231314@gmail.com
% Joint China-Cuba LAB, UESTC, Hangzhou Dianzi University
% License: GNU General Public License v3.0 (see LICENSE file)

function [x,lb,ub]=recover_fixed_para(xsub,lbsub,ubsub,parainfo)

% if ~isempty(xsub)
x=parainfo.x;
x(~parainfo.xIndices.fixed)=xsub;
% else
%     x=[];
% end


if ~isempty(lbsub)
    lb=parainfo.lb;
    lb(~parainfo.xIndices.fixed)=lbsub;
else
    lb=[];
end
if ~isempty(ubsub)
    ub=parainfo.ub;
    ub(~parainfo.xIndices.fixed)=ubsub;
else
    ub=[];
end

end



%remove_fixed_para