% Author: Ying Wang, Min Li
% Create Time: 2025
% Copyright (c): 2020-2025 Ying Wang, yingwangrigel@gmail.com,
%                Min Li, minli.231314@gmail.com
% Joint China-Cuba LAB, UESTC, Hangzhou Dianzi University
% License: GNU General Public License v3.0 (see LICENSE file)

function shat = shat_component(self, component)

switch component
    case {'xi'}
        [shat] = pred_s_xi(self.f, self.para_xi);
    case {'alpha'}
        [shat] = pred_s_alpha(self.f, self.ks, self.kh, self.para_alpha); %,comp gaussian_kernel(f,mu,sigma).*para_alpha.kernel.alpha.';

    otherwise
        error('the component is not exsit, only xi and alpha');
end

end

% para_bs.betaxy(2:end)=0;
% para_bs.gammaxy(2:end)=0;

