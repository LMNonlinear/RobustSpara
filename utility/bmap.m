% Author: Ying Wang, Min Li
% Create Time: 2025
% Copyright (c): 2020-2025 Ying Wang, yingwangrigel@gmail.com,
%                Min Li, minli.231314@gmail.com
% Joint China-Cuba LAB, UESTC, Hangzhou Dianzi University
% License: GNU General Public License v3.0 (see LICENSE file)

function A = bmap(mA, idx, A)
% compare object type of mA and A, if class(A(i)) is not the same as class(mA(i)) then cast it
if ~strcmpi(class(mA), class(A))
    A=create_empty_objarray(size(A),mA);
    % a(2)=mA(1);
    % A = a;
    % cast(A, class(mA))
    % A=(class(mA)).empty(size(A));
end

A(idx) = mA;
end

% function A=create_empty_array_obj(sz,a)
% 
% % A(2)=a(1);
% % A=repmat(A(1),sz);
% A=repmat(get_empty_obj(a),sz);
% end


% function A=get_empty_obj(a)
% if length(a)>1
%     A(2)=a(1);
% else
%     A(2)=a;
% end
% end



