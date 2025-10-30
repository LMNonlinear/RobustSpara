% Author: Ying Wang, Min Li
% Create Time: 2025
% Copyright (c): 2020-2025 Ying Wang, yingwangrigel@gmail.com,
%                Min Li, minli.231314@gmail.com
% Joint China-Cuba LAB, UESTC, Hangzhou Dianzi University
% License: GNU General Public License v3.0 (see LICENSE file)

function [bshat] = pred_bs(f, para_bs,fxfy)

para_bs.kernel.para.h(isnan(para_bs.kernel.para.h)) = 0;
bshat=zeros(size(f,1),size(f,1));
switch para_bs.type

    case 'tstudent'
        % bshat=para_bs.kernel.eval(fxfy);
        % bshat=atriu2full(bshat,true);
        % bshat=tril2full(bshat,false,false);
        [bshat(para_bs.idx_fxfy)]=para_bs.kernel.eval(fxfy);
        bshat=tril2full(bshat,false,false);
    otherwise
        error('no such method')
        
end



end
