% Author: Ying Wang, Min Li
% Create Time: 2025
% Copyright (c): 2020-2025 Ying Wang, yingwangrigel@gmail.com,
%                Min Li, minli.231314@gmail.com
% Joint China-Cuba LAB, UESTC, Hangzhou Dianzi University
% License: GNU General Public License v3.0 (see LICENSE file)

function [index,res]=find_closest(A,a,n)
if nargin<3
    n=1;
end
% A=[34.8 31 29 26.7 39.5];%dummy data
% a=33;
% [index,res]=find_closest(A,a)
if n==1
    index=arrayfun(@(x) minidx(abs(A-x)),a);
    res=A(index)-a;
else
    index=arrayfun(@(x) sortidx(abs(A-x),n),a,'UniformOutput',false);
    index=cell2mat(index);
    res=A(index)-a;
end


    function idx=minidx(x)
        [~,idx]=min(x);
    end
    function idx=sortidx(x,n)
        [~,idx]=sort(x);
        idx=idx(1:n);
    end

end