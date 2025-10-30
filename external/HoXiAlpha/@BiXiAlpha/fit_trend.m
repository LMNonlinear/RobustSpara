% Author: Ying Wang, Min Li
% Create Time: 2025
% Copyright (c): 2020-2025 Ying Wang, yingwangrigel@gmail.com,
%                Min Li, minli.231314@gmail.com
% Joint China-Cuba LAB, UESTC, Hangzhou Dianzi University
% License: GNU General Public License v3.0 (see LICENSE file)

function trend = fit_trend(self)
% fit trend later use this trend to get the parameter of xi

%% only fit xi ignore the peak, nonpara
if strcmpi(self.method_xi, 'robust') % % robust detrend
    % [~,~,trend] = RobustDetrend(log10(self.s),10,0.975,self.f);%LT= trenddecomp(self.s);
    [~, ~, trend] = RobustDetrend(convert_to_double(mean(log(self.s), 2)), 3, 0.975, self.f); %LT= trenddecomp(self.s);

elseif strcmpi(self.method_xi, 'linear') % % linear detrend
    n = 1;
    trend = detrend(log(mean(self.s, 2)), n);
end

% fit trend in log scale, addtive model in the raw scale
trend = exp(trend(:));
% figure(101);plot(log10([self.s,trend]));hold on

end

