% Author: Ying Wang, Min Li
% Create Time: 2025
% Copyright (c): 2020-2025 Ying Wang, yingwangrigel@gmail.com,
%                Min Li, minli.231314@gmail.com
% Joint China-Cuba LAB, UESTC, Hangzhou Dianzi University
% License: GNU General Public License v3.0 (see LICENSE file)

function rsquared = calc_goodnessOfFit(y, yfit)
y=mean(y,2,'omitnan');


% Calculate the goodness of fit
SSE = sum((abs(y - yfit)).^2,"all",'omitnan');           % Sum of Squared Errors
SST = sum((abs(y - mean(y,1,'omitnan'))).^2,"all",'omitnan');        % Total Sum of Squares
rsquared = 1 - SSE/SST;             % R-squared

end











% % Calculate the goodness of fit
% SSE = sum((abs(y - yfit)).^2,"all",'omitnan');           % Sum of Squared Errors
% SST = sum((abs(y - mean(y,1,'omitnan'))).^2,"all",'omitnan');        % Total Sum of Squares
% rsquared = 1 - SSE/SST;             % R-squared
