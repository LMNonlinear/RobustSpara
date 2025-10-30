% Author: Ying Wang, Min Li
% Create Time: 2025
% Copyright (c): 2020-2025 Ying Wang, yingwangrigel@gmail.com,
%                Min Li, minli.231314@gmail.com
% Joint China-Cuba LAB, UESTC, Hangzhou Dianzi University
% License: GNU General Public License v3.0 (see LICENSE file)

function surf(self)
    bhat = pred_bs(self);
    bsre = real(bhat);
    bsim = imag(bhat);

    subplot(2, 1, 1)
    p1 = surf(self.fx, self.fy, abs(self.bs)); hold on;
    view(0, 90)
    set(p1, 'edgecolor', 'none')
    shading interp
    grid off
    colormap(flip(hot))
    colorbar

    subplot(2, 1, 2)
    p2 = surf(self.fx, self.fy, sqrt(bsre .^ 2 + bsim .^ 2));
    view(0, 90)
    set(p2, 'edgecolor', 'none')
    shading interp
    grid off
    colormap(flip(hot))
    colorbar
end

