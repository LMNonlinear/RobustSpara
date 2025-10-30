% Author: Ying Wang, Min Li
% Create Time: 2025
% Copyright (c): 2020-2025 Ying Wang, yingwangrigel@gmail.com,
%                Min Li, minli.231314@gmail.com
% Joint China-Cuba LAB, UESTC, Hangzhou Dianzi University
% License: GNU General Public License v3.0 (see LICENSE file)


function fit(self)
self.bs = get_bs2fit(self);
self.init_para();

switch lower(self.optimization)
    case 'gmm'
        self.GMM();
    case 'nonlinearoptimization'
        self.NonlinearOptimization();
    otherwise
        error('no such method')
end
disp(['finish fit for name: [',num2str(self.name),']']);
% disp(repmat(newline, 1, 1));

end
