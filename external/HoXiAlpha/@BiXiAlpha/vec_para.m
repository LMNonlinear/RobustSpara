% Author: Ying Wang, Min Li
% Create Time: 2025
% Copyright (c): 2020-2025 Ying Wang, yingwangrigel@gmail.com,
%                Min Li, minli.231314@gmail.com
% Joint China-Cuba LAB, UESTC, Hangzhou Dianzi University
% License: GNU General Public License v3.0 (see LICENSE file)

function x=vec_para(self,part,ispad2max)
if nargin<3 ||isempty(ispad2max)
    ispad2max=true;
end

switch part
    case {'k'}
        x=[self.ks;self.kh];
    case {'xi'}
        % switch self.para_xi.type
        %     case{'tstudent'}
        %         x=self.para_xi.kernel.vpara();
        %     otherwise
        %         error('no such method')
        % end
        x=self.para_xi.kernel.vpara();
    case{'alpha'}
        x=self.para_alpha.kernel.vpara();
    case{'bs'}
        x=self.para_bs.kernel.vpara();
    otherwise
        error('no such part')
end
end
