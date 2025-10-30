% Author: Ying Wang, Min Li
% Create Time: 2025
% Copyright (c): 2020-2025 Ying Wang, yingwangrigel@gmail.com,
%                Min Li, minli.231314@gmail.com
% Joint China-Cuba LAB, UESTC, Hangzhou Dianzi University
% License: GNU General Public License v3.0 (see LICENSE file)

function S=set_field(S,field,val,index)
if ischar(field)||isstring(field)
    field={field};
end

if nargin<4
    if ~isempty(S)
        index=1:numel(S);
    else
        index=1;
    end
end
for ifield=1:numel(field)
    for iidx=1:numel(index)
        if numel(val)>1 &&  ~ischar(val) && ~isstring(val) &&numel(val)==numel(index)
            if ismatrix(val)
                S = setfield(S,{index(iidx)},field{ifield},val(index(iidx)));
            elseif iscell(val)
                S = setfield(S,{index(iidx)},field{ifield},val{index(iidx)});
            end
        else
            S = setfield(S,{index(iidx)},field{ifield},val);
        end
    end
end
end