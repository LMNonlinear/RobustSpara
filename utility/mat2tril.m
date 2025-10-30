% Author: Ying Wang, Min Li
% Create Time: 2025
% Copyright (c): 2020-2025 Ying Wang, yingwangrigel@gmail.com,
%                Min Li, minli.231314@gmail.com
% Joint China-Cuba LAB, UESTC, Hangzhou Dianzi University
% License: GNU General Public License v3.0 (see LICENSE file)

function [mattril,idx]=mat2tril(data,k,isvec,onlyidx,isfirstslice)
% examp: 10*10*100, get lower triangle
% Yingwang 5/10/2020
%  ind= 0 is the main diagonal
%  ind > 0 is above the main diagonal
%  ind < 0 is below the main diagonal.
if nargin<5 || isempty(isfirstslice)
    isfirstslice=false;
end
if nargin<4 || isempty(onlyidx)
    onlyidx=false;
end
if nargin<3 || isempty(isvec)
    isvec=true;
end
if nargin<2 || isempty(k)
    k=0;
end
[r,c,p]=size(data);

if isvec
    % idx=tril(true(r,c),k);
    % idx=find(idx);
    % idx=idx(:)+r*c*(0:p-1);
    % if ~onlyidx
    %     mattril=reshape(data(idx),[],p);
    % end
    if ~isfirstslice
        idx=tril(true(r,c),k);
        idx=find(idx);
        idx=idx(:)+r*c*(0:p-1);
        if ~onlyidx
            mattril=reshape(data(idx),[],p);
        end
    else
        idx=tril(true(r,c),k);
        idx=find(idx);
        if ~onlyidx
            data=reshape(data,[r*c,p]);
            mattril=reshape(data(idx,:),[],p);
        end
    end
else
    if ~onlyidx
        mattril=zeros(size(data));
        for i=1:p
            mattril(:,:,i)=tril(data(:,:,i),k);
        end
    end
    if nargout>1
        idx=tril(true(r,c),k);
        idx=repmat(idx,[1,1,p]);
    end
end
if onlyidx
    mattril=[];
end



end














