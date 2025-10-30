% Author: Ying Wang, Min Li
% Create Time: 2025
% Copyright (c): 2020-2025 Ying Wang, yingwangrigel@gmail.com,
%                Min Li, minli.231314@gmail.com
% Joint China-Cuba LAB, UESTC, Hangzhou Dianzi University
% License: GNU General Public License v3.0 (see LICENSE file)

function shat=get_shat_component(self,component)
    switch component
        case {'xi'}
            shat=pred_s_xi(self.f,self.para_xi);
        case {'alpha'}
            shat=pred_s_alpha(self.f,self.kh,self.para_alpha);
        otherwise
            error('the component is not exsit, only xi and alpha');
    end

end
