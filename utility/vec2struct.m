% Author: Ying Wang, Min Li
% Create Time: 2025
% Copyright (c): 2020-2025 Ying Wang, yingwangrigel@gmail.com,
%                Min Li, minli.231314@gmail.com
% Joint China-Cuba LAB, UESTC, Hangzhou Dianzi University
% License: GNU General Public License v3.0 (see LICENSE file)

function s = vec2struct(vec, info)
% Recover a structure from a vector.
s = struct();
start_idx = 1;

if isvector(vec)

    for i = 1:size(info.fields, 1)
        field = info.fields{i, 1};
        sz = info.fields{i, 2};
        num_elements = prod(sz);
        s.(field) = reshape(vec(start_idx:start_idx + num_elements - 1), sz);
        start_idx = start_idx + num_elements;
    end

else

    for i = 1:size(info.fields, 1)
        field = info.fields{i, 1};
        sz = info.fields{i, 2};
        num_elements = prod(sz);

        if start_idx + num_elements - start_idx .* size(vec, 2) > prod(sz)
            sz = [sz(:); size(vec, 2)]';
        end

        s.(field) = reshape(vec(start_idx:start_idx + num_elements - 1, :), sz);
        start_idx = start_idx + num_elements;
    end

end

end
