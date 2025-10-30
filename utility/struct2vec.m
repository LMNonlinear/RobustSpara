% Author: Ying Wang, Min Li
% Create Time: 2025
% Copyright (c): 2020-2025 Ying Wang, yingwangrigel@gmail.com,
%                Min Li, minli.231314@gmail.com
% Joint China-Cuba LAB, UESTC, Hangzhou Dianzi University
% License: GNU General Public License v3.0 (see LICENSE file)

function [vec, info] = struct2vec(s)
% Convert a structure to a vector.
fields = fieldnames(s);
info.fields = cell(length(fields), 2);
info.vecLength = 0;

% First loop to gather information and calculate total length
for i = 1:length(fields)
    field = s.(fields{i});
    info.fields{i, 1} = fields{i};
    info.fields{i, 2} = size(field);
    info.vecLength = info.vecLength + numel(field);
end

% Preallocate vec
vec = zeros(info.vecLength, 1);

% Second loop to fill vec
index = 1;
for i = 1:length(fields)
    field = s.(fields{i});
    numElems = numel(field);
    vec(index : index + numElems - 1) = field(:);
    index = index + numElems;
end
end

% function [vec, info] = struct2vec(s)
% % Convert a structure to a vector.
% fields = fieldnames(s);
% info = cell(length(fields), 2);
% vec = [];
% for i = 1:length(fields)
%     field = s.(fields{i});
%     info{i, 1} = fields{i};
%     info{i, 2} = size(field);
%     vec = [vec; field(:)];  % Concatenate into a single vector
% end
% end

