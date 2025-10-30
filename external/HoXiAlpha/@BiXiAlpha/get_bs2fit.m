% Author: Ying Wang, Min Li
% Create Time: 2025
% Copyright (c): 2020-2025 Ying Wang, yingwangrigel@gmail.com,
%                Min Li, minli.231314@gmail.com
% Joint China-Cuba LAB, UESTC, Hangzhou Dianzi University
% License: GNU General Public License v3.0 (see LICENSE file)


% recover compacted data back to full matrix
function bs = get_bs2fit(self)

switch self.compact
    case 'wedge'
        % bs=mat2wedge(self.bs,false);
        bs = wedge2full(self.bs, true);
    case 'tril'
        % bs=mat2atriu(self.bs,0,false);
        % bs=mat2tril(self.bs,0,false);
        bs = tril2full(self.bs, 0, false);
    case 'quad1'
        bs = self.bs;
    case 'full'
        bs = self.bs;
end

end

