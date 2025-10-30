% Author: Ying Wang, Min Li
% Create Time: 2025
% Copyright (c): 2020-2025 Ying Wang, yingwangrigel@gmail.com,
%                Min Li, minli.231314@gmail.com
% Joint China-Cuba LAB, UESTC, Hangzhou Dianzi University
% License: GNU General Public License v3.0 (see LICENSE file)

function plot(self)

    sest=self.shat;
    subplot(2,1,1)
    plot(self.f,self.s);hold on

    subplot(2,1,2)
    plot(self.f,sum(sest,2));hold on
    % plot(self.f,sest);
end
