% Author: Ying Wang, Min Li
% Create Time: 2025
% Copyright (c): 2020-2025 Ying Wang, yingwangrigel@gmail.com,
%                Min Li, minli.231314@gmail.com
% Joint China-Cuba LAB, UESTC, Hangzhou Dianzi University
% License: GNU General Public License v3.0 (see LICENSE file)


function bs=compact_bs(bs,compact)
if ischar(compact) || isstring(compact)
    switch compact
        case 'wedge'
            bs=mat2wedge(bs,false);
        case 'tril'
            % bs=mat2atriu(self.bs,0,false);
            bs=mat2tril(bs,0,false);
        case 'quad1'
            bs=bs;
        case 'full'
            bs=bs;
        otherwise
            error('no such compact')
    end
else
    bs=bs(compact);
end

end