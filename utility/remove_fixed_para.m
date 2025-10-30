% Author: Ying Wang, Min Li
% Create Time: 2025
% Copyright (c): 2020-2025 Ying Wang, yingwangrigel@gmail.com,
%                Min Li, minli.231314@gmail.com
% Joint China-Cuba LAB, UESTC, Hangzhou Dianzi University
% License: GNU General Public License v3.0 (see LICENSE file)

function [x,lb,ub,parainfo]=remove_fixed_para(x,lb,ub)
nvar=length(x);
if nvar~=length(lb)||nvar~=length(ub)
    error('x, lb, and ub must have the same length');
end
parainfo.xIndices = classifyBoundsOnVars(lb,ub,nvar,true);
parainfo.x=x;
parainfo.lb=lb;
parainfo.ub=ub;

x(parainfo.xIndices.fixed)=lb(parainfo.xIndices.fixed);

x=x(~parainfo.xIndices.fixed);
lb=lb(~parainfo.xIndices.fixed);
ub=ub(~parainfo.xIndices.fixed);

end


%recover_fixed_para