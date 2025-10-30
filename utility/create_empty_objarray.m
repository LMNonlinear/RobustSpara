% Author: Ying Wang, Min Li
% Create Time: 2025
% Copyright (c): 2020-2025 Ying Wang, yingwangrigel@gmail.com,
%                Min Li, minli.231314@gmail.com
% Joint China-Cuba LAB, UESTC, Hangzhou Dianzi University
% License: GNU General Public License v3.0 (see LICENSE file)

function A = create_empty_objarray(sz, a)

% A(2)=a(1);
% A=repmat(A(1),sz);
A = repmat(get_empty_obj(a), sz);
end
