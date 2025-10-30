% Author: Ying Wang, Min Li
% Create Time: 2025
% Copyright (c): 2020-2025 Ying Wang, yingwangrigel@gmail.com,
%                Min Li, minli.231314@gmail.com
% Joint China-Cuba LAB, UESTC, Hangzhou Dianzi University
% License: GNU General Public License v3.0 (see LICENSE file)

function init_para_alpha(self)

% para_alpha=self.para_alpha;
% if strcmpi(self.para_alpha.type,'gaussian')
%     para_alpha.kernel.para.mu=11;% only for xi and alpha(harmonic)
%     para_alpha.kernel.para.sigma=1;% only for xi and alpha
%     para_alpha.kernel.alpha=[ones(1,self.ks+self.kh)];% alpha might have relation between harmonic, add later
%     % now only for spontaneous alpha, intermodulation subharmonic is not
%     % considered yet, can be add easily
%     para_alpha.lb_mu=4;
%     para_alpha.ub_mu=30;
%     % para_alpha.lb_mu=5;
%     % para_alpha.ub_mu=15;
% else
%     error('no such method')
% end

switch self.para_alpha.type
    case 'tstudent'

        if strcmpi(self.para_alpha.type, 'tstudent') && strcmpi(self.peak_relation, 'harmonic')

            if ~isfield(self.para_alpha, 'kernel') || isempty(self.para_alpha.kernel) || ~strcmpi(self.para_alpha.kernel.type, self.kernel_type) || isempty(self.para_alpha.kernel.para) || length(self.para_alpha.kernel.para.h) ~= self.kmax
                self.para_alpha.kernel = Kernel(self.kernel_type);
                self.para_alpha.kernel.info.no_zero_peak = self.para_alpha.no_zero_peak;
                self.k0 = double(~self.para_alpha.kernel.info.no_zero_peak);

                self.para_alpha.kernel.para.h = [ones(1, self.kmax)];
                self.para_alpha.kernel.para.mu = 11;
                self.para_alpha.kernel.para.sigma = 1.39; %10 %
                self.para_alpha.kernel.para.nu = 2; % 2 %60
                self.para_alpha.kernel.para.d = 6.14; %3
                self.para_alpha.kernel.para.b = 1;
                self.para_alpha.kernel.info.ks = self.ks;
                self.para_alpha.kernel.info.kh = self.kh;
                self.para_alpha.kernel.info.idx_sigma = [];
                self.para_alpha.kernel.info.df = self.df;
                self.para_alpha.kernel.info.taperinfo = self.taperinfo;
            end

        elseif strcmpi(self.para_alpha.type, 'tstudent') && strcmpi(self.peak_relation, 'free')

            if ~isfield(self.para_alpha, 'kernel') || isempty(self.para_alpha.kernel) || ~strcmpi(self.para_alpha.kernel.type, self.kernel_type) || isempty(self.para_alpha.kernel.para) || length(self.para_alpha.kernel.para.h) ~= self.kmax
                self.para_alpha.kernel = Kernel(self.kernel_type);
                self.para_alpha.kernel.info.no_zero_peak = self.para_alpha.no_zero_peak;
                self.k0 = double(~self.para_alpha.kernel.info.no_zero_peak);

                self.para_alpha.kernel.para.h = [ones(1, self.kmax)];
                self.para_alpha.kernel.para.mu = 11 .* ones(self.kmax, 1);
                self.para_alpha.mu_init = self.para_alpha.kernel.para.mu;
                self.para_alpha.kernel.para.sigma = 1.39 .* ones(self.kmax, 1); %10 %
                self.para_alpha.kernel.para.nu = 2 .* ones(self.kmax, 1); % 2 %60
                self.para_alpha.kernel.para.d = 6.14 .* ones(self.kmax, 1); %3
                self.para_alpha.kernel.para.b = 1 .* ones(self.kmax, 1);
                self.para_alpha.kernel.info.ks = self.ks;
                self.para_alpha.kernel.info.kh = self.kh;
                self.para_alpha.kernel.info.idx_sigma = [];
                self.para_alpha.kernel.info.df = self.df;
                self.para_alpha.kernel.info.taperinfo = self.taperinfo;
            end

        else
            error('no such peak_relation')
        end

    case 'gaussian'

        if ~isfield(self.para_alpha, 'kernel') || isempty(self.para_alpha.kernel) || ~strcmpi(self.para_alpha.kernel.type, self.kernel_type) || isempty(self.para_alpha.kernel.para) || length(self.para_alpha.kernel.para.h) ~= self.kmax
            self.para_alpha.kernel = Kernel('gaussian_kernel');
            self.para_alpha.kernel.info.no_zero_peak = self.para_alpha.no_zero_peak;
            self.k0 = double(~self.para_alpha.kernel.info.no_zero_peak);

            self.para_alpha.kernel.para.h = [ones(1, self.kmax)];
            self.para_alpha.kernel.para.mu = 11 .* ones(self.kmax, 1);
            self.para_alpha.mu_init = self.para_alpha.kernel.para.mu;
            self.para_alpha.kernel.para.sigma = 0.4 .* ones(self.kmax, 1); %10 %
            self.para_alpha.kernel.info.ks = self.ks;
            self.para_alpha.kernel.info.kh = self.kh;
            self.para_alpha.kernel.info.idx_sigma = [];
            self.para_alpha.kernel.info.df = self.df;
            self.para_alpha.kernel.info.taperinfo = self.taperinfo;
        end

    otherwise
        error('no such method')
end

end

