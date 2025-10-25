classdef Kernel %< handle & matlab.mixin.Copyable

    properties
        kfun
        type
        para % a struct of parameters
        lb
        ub
        fix
        t
        % res
        % jac
        info
        use_info=false;
        hist
        % D

    end

    properties (Hidden)
        vpara % a vector of parameters
        vlb % a vector of lower bounds of parameters
        vub
        vfix
        vt
        % vres
        % vjac
        %store the information of parameters structure
        para_info
        lb_info
        ub_info
        fix_info
        t_info
        % res_info
        % jac_info
    end

    methods

        function self = Kernel(type, para)

            if nargin == 0
                self.type = 'gaussian_kernel';
                self.para = struct('h', 1, 'mu', 0, 'sigma', 1);
            elseif nargin == 1
                self.type = type;
                self = self.init_para(type);
            elseif nargin == 2
                self.type = type;
                self.para = para;
            end

            % self.hist.update_dist=false;
        end

        %% set methods
        function self = set.para(self, para)
            self.para = para;
            [~, self.para_info] = struct2vec(para);
        end

        function self = set.lb(self, lb)
            self.lb = lb;
            [~, self.lb_info] = struct2vec(lb);
        end

        function self = set.ub(self, ub)
            self.ub = ub;
            [~, self.ub_info] = struct2vec(ub);
        end

        function self = set.t(self, t)
            self.t = t;
            [~, self.t_info] = struct2vec(t);
        end

        function self = set.vpara(self, value)
            self.vpara = value;
            % set all parameters from a vector
            self.para = vec2struct(value, self.para_info);
        end

        function self = set.vlb(self, value)
            self.vlb = value;

            if isempty(self.lb_info)
                self.lb_info = self.para_info;
            end

            % set all lower bounds of parameters from a vector
            self.lb = vec2struct(value, self.lb_info);
        end

        function self = set.vub(self, value)
            self.vub = value;

            if isempty(self.ub_info)
                self.ub_info = self.para_info;
            end

            % set all upper bounds of parameters from a vector
            self.ub = vec2struct(value, self.ub_info);
        end

        function self = set.vfix(self, value)
            self.vfix = value;

            if isempty(self.fix_info)
                self.fix_info = self.para_info;
            end

            % set all fix value of parameters from a vector
            self.fix = vec2struct(value, self.fix_info);
        end

        function self = set.vt(self, value)
            self.vt = value;

            if isempty(self.t_info)
                self.t_info = self.para_info;
            end

            % set all t value of parameters from a vector
            self.t = vec2struct(value, self.t_info);
        end

        % function self = set.vres(self, value)
        %     self.vres = value;

        %     if isempty(self.res_info)
        %         self.res_info = self.para_info;
        %     end

        %     % set all res value of parameters from a vector
        %     self.res = vec2struct(value, self.res_info);
        % end

        % function self = set.vjac(self, value)
        %     self.vjac = value;

        %     if isempty(self.jac_info)
        %         self.jac_info = self.para_info;
        %     end

        %     % set all jac value of parameters from a vector
        %     self.jac = vec2struct(value, self.jac_info);
        % end

        function self = set.type(self, value)
            self.type = value;
            self.kfun = str2func(['Kernel.', value]); %#ok<*MCSUP>
        end

        %% get methods
        function vpara = get.vpara(self)
            % concatenate all parameters to a vector
            [vpara] = struct2vec(self.para); %,self.para_info
        end

        function vlb = get.vlb(self)
            % concatenate all lower bounds of parameters to a vector
            [vlb] = struct2vec(self.lb); %,self.lb_info
        end

        function vub = get.vub(self)
            % concatenate all upper bounds of parameters to a vector
            [vub] = struct2vec(self.ub); %,self.ub_info
        end

        function vfix = get.vfix(self)

            if ~isempty(self.fix)
                % concatenate all fix value of parameters to a vector
                [vfix] = struct2vec(self.fix); %,self.fix_info
            end

        end

        function vt = get.vt(self)
            % concatenate all t value of parameters to a vector
            [vt] = struct2vec(self.t); %,self.t_info
        end

        % function vres = get.vres(self)
        %     % concatenate all res value of parameters to a vector
        %     [vres] = struct2vec(self.res); %,self.res_info
        % end

        % function vjac = get.vjac(self)
        %     % concatenate all jac value of parameters to a vector
        %     warning('vjac is not used, direct vectorized jacobian is no use')
        %     [vjac] = struct2vec(self.jac); %,self.jac_info
        % end

        %% other methods
        function [y, ycomp] = eval(self, x)
            % tic
            % feed structure para as cell format to kfun to get y
            para = struct2cell(self.para); %#ok<*PROPLC>
            if self.use_info||strcmpi(self.type,'harmonic_t_flattop')
                info = struct2cell(self.info);
            end
            if nargout == 1
                if self.use_info||strcmpi(self.type,'harmonic_t_flattop')
                    y = self.kfun(x, para{:}, info{:});
                else
                    y = self.kfun(x, para{:});
                end
            else
                if self.use_info ||strcmpi(self.type,'harmonic_t_flattop')
                    [y, ycomp] = self.kfun(x, para{:}, info{:});
                else
                    [y, ycomp] = self.kfun(x, para{:});
                end
            end

            % toc
        end

        function T = check(self)
            T = [self.lb, self.para, self.ub];
            % check if the parameters are within the bounds
            % 0: not within bounds, 1: within bounds

            % T=[self.vlb,self.vpara,self.vub];
            if ~isempty(self.vlb) && ~isempty(self.vpara) && ~isempty(self.vub)
                paras = [self.vlb, self.vpara, self.vub];
                tf = [paras(:, 1) <= paras(:, 2), paras(:, 2) <= paras(:, 3)];

                if any(tf(:) == 0)
                    disp('out of bounds')
                else
                    disp('within bounds')
                end

            end

        end

        function self = init_para(self, type)

            switch type
                case 'gaussian_kernel'
                    self.para = struct('h', 1, 'mu', 0, 'sigma', 1);
                case 'lorentzian_kernel'
                    self.para = struct('h', 1, 'mu', 0,'k',1, 'chi', 1);
                case 'linear_kernel'
                    self.para = struct('h', 1);
                case 'polynomial_kernel'
                    self.para = struct('h', 1, 'd', 2);
                case 'exponential_kernel'
                    self.para = struct('h', 1, 'beta', 1);
                    % case 'power_exponential_kernel'
                    %     self.para=struct('h',1,'beta',1,'d',2);
                case 'studentt_kernel'
                    self.para = struct('A', 1, 'B', 1, 'nu', 1, 'U', 0);
                    % self.para=struct('A',600,'B',9,'nu',3.2,'U',10);
                case 'laplacian_kernel'
                    self.para = struct('h', 1, 'beta', 1);
                case 'periodic_kernel'
                    self.para = struct('h', 1, 'beta', 1, 'period', 1);
                case 'matern_kernel'
                    self.para = struct('h', 1, 'beta', 1, 'nu', 1);
                case 'modified_gaussian_kernel'
                    self.para = struct('h', 1, 'mu', 0, 'sigma', 1, 'R', 10);
                case 'bs_kernel'
                    self.para = struct('h', 1, 'mu', [0, 0; 10, 10], 'sigma', [1, 1, -0.2]);
                case 'bs_kernel_flat'
                    self.para = struct('h', 1, 'mu', [0, 0; 10, 10], 'sigma', [1, 1, -0.2], 'R', 10);
                case 'bs_kernel_t'
                    self.para = struct('h', 1, 'mu', [0, 0; 10, 10], 'sigma', [1, 1, -0.2], 'nu', 10);
                case 'studentt_waveform'
                    self.para = struct('h', 1, 'mu', [0, 0; 10, 10], 'sigma', [2, 2, -0.2], 'nu', 10);
                case 'harmonic_t'
                    self.para = struct('h', [1, 0; 0, 2], 'mu', [0; 10], 'sigma', [10, 1; 1, 1], 'nu', [1; 1]);
                    self.info = struct('ks', 0, 'kh', 1, 'idx_sigma', [1, 0; 0, 2], 'no_zero_peak', false, 'df', 1);
                case 'harmonic_t_flattop'
                    self.use_info=true;
                    self.para = struct('h', [1, 0; 0, 2], 'mu', [0; 10], 'sigma', [10, 1; 1, 1], 'nu', [1; 1], 'd', [1, 1], 'b', [1, 1]);
                    self.info = struct('ks', 0, 'kh', 1, 'idx_sigma', [1, 0; 0, 2]);
                case 'harmonic_t_free'
                    self.use_info=true;
                    self.para = struct('h', [1, 0; 0, 2], 'mu', [0; 10], 'sigma', [10, 1; 1, 1], 'nu', [1; 1], 'd', [1, 1], 'b', [1, 1]);
                    self.info = struct('ks', 0, 'kh', 1, 'idx_sigma', [1, 0; 0, 2]);
                case 'harmonic_t_tapers'
                    self.para = struct('h', [1, 0; 0, 2], 'mu', [0; 10], 'sigma', [10, 1; 1, 1], 'nu', [1; 1]);
                    self.info = struct('ks', 0, 'kh', 1, 'idx_sigma', [1, 0; 0, 2]);
                otherwise
                    error('no such method')
            end

        end

    end

    methods (Static)

        function para = setVectorField(para, value)
            % set all parameters of fields from a vector
            % 1. get the field names
            fnames = fieldnames(para);
            % 2. get the length of each field
            flength = cellfun(@(x) length(para.(x)), fnames);
            % 3. set the value of each field from the vector with vectorization
            for i = 1:length(flength)

                if i == 1
                    para.(fnames{i}) = value(1:flength(i));
                else
                    para.(fnames{i}) = value(sum(flength(1:i - 1)) + 1:sum(flength(1:i)));
                end

            end

        end

        function [y,ycomp] = gaussian_kernel(x, h, mu, sigma)
            if nargout==1
                y = sum(h' .* exp(- (x - mu') .^ 2 ./ (2 * sigma' .^ 2)),2);
            elseif nargout==2
                ycomp = h' .* exp(- (x - mu') .^ 2 ./ (2 * sigma' .^ 2));
                y=sum(ycomp,2);
            end
        end

        function [y,ycomp] = lorentzian_kernel(x, h, mu, k, chi)
            if nargout==1
                y = sum(h' ./ (k + (x-mu') .^ (chi')),2);
            elseif nargout==2
                ycomp = h' ./ (k + (x-mu') .^ (chi'));
                y=sum(ycomp,2);
            end
        end

        function y = linear_kernel(x, h)
            y = h * x;
        end

        function y = polynomial_kernel(x, h, d)
            y = h * x .^ d;
        end

        function y = exponential_kernel(x, h, beta)
            y = h * exp(-beta * x);
        end

        % function y=power_exponential_kernel(x,h,beta,nu)
        %     y=h*;
        % end
        function y = studentt_kernel(x, A, B, nu, U)
            % follow the definition of student-t kernel in the Pascual's paper
            y = A ./ ((1.0 + ((x - U) ./ B) .^ 2) .^ nu);
        end

        function y = laplacian_kernel(x, h, beta)
            y = h * exp(-beta * abs(x));
        end

        function y = periodic_kernel(x, h, beta, period)
            y = h * exp(-2 * beta * sin(pi * abs(x) / period) .^ 2);
        end

        function y = matern_kernel(x, h, beta, nu)
            y = h * (1 + sqrt(2 * nu) / beta * abs(x)) .* exp(-sqrt(2 * nu) / beta * abs(x));
        end

        function y = modified_gaussian_kernel(x, h, mu, sigma, R)
            y = h * (1 + R) ./ (R + exp((x - mu) .^ 2 / (2 * sigma ^ 2)));
        end

        function y = bs_kernel(x, h, mu, sigma)
            % x is 2d meshgrid , h is the height of kernels, mu is the center of kernels, sigma is the width of kernels, all kernel share the same sigma
            % h is Nk*1 vector, mu is Nk*2 matrix ,sigma is 3*1 matrix, x is Nx*2 matrix, y is Nx*1 vector, Nk is the number of kernels, Nx is the number of data points

            % f(x,y)=h*exp(-1/2*(x-mu)'*Sigma^(-1)*(x-mu))
            % recover covariance matrix from sigma to Sigma, Sigma is 2*2 matrix
            Sigma = [sigma(1), conj(sigma(3)); sigma(3), sigma(2)];
            Nx = size(x, 1);
            x = reshape(x, [], 2);

            % inverse of Sigma
            % invSigma = inv(Sigma);
            % calculate Mahalanobis distance
            D = pdist2(x, mu, 'mahalanobis', Sigma);

            % calculate kernel for each data point and each kernel
            y = sum(h .* exp(-0.5 * D .^ 2), 2);
            y = reshape(y, Nx, Nx);

        end

        function y = bs_kernel_flat(x, h, mu, sigma, R)
            % x is 2d meshgrid , h is the height of kernels, mu is the center of kernels, sigma is the width of kernels, all kernel share the same sigma
            % h is Nk*1 vector, mu is Nk*2 matrix ,sigma is 3*1 matrix, x is Nx*2 matrix, y is Nx*1 vector, Nk is the number of kernels, Nx is the number of data points

            % f(x,y)=h*exp(-1/2*(x-mu)'*Sigma^(-1)*(x-mu))
            % recover covariance matrix from sigma to Sigma, Sigma is 2*2 matrix
            Sigma = [sigma(1), conj(sigma(3)); sigma(3), sigma(2)];
            Nx = size(x, 1);
            x = reshape(x, [], 2);

            % inverse of Sigma
            % invSigma = inv(Sigma);
            % calculate Mahalanobis distance
            D = pdist2(x, mu, 'mahalanobis', Sigma);

            % calculate kernel for each data point and each kernel
            y = sum(h .* (1 + R) ./ (R + exp(0.5 * D .^ 2)), 2);
            y = reshape(y, Nx, Nx);

        end

        function y = bs_kernel_t(x, h, mu, sigma, nu)
            % x is 2d meshgrid, h is the height of kernels, mu is the center of kernels, sigma is the width of kernels, all kernel share the same sigma, nu is the degree of freedom
            % h is Nk*1 vector, mu is Nk*2 matrix ,sigma is 3*1 matrix, nu is scalar, x is Nx*2 matrix, y is Nx*1 vector, Nk is the number of kernels, Nx is the number of data points

            % f(x,y)=h*(1+D^2/nu)^(-(nu+d)/2)
            % recover covariance matrix from sigma to Sigma, Sigma is 2*2 matrix
            Sigma = [sigma(1), conj(sigma(3)); sigma(3), sigma(2)];

            if ndims(x) == 3
                Nx = size(x, 1);
                x = reshape(x, [], 2);
            end

            % calculate Mahalanobis distance
            D = pdist2(x, mu, 'mahalanobis', Sigma);

            % calculate kernel for each data point and each kernel
            y = sum(h .* (1 + D .^ 2 / nu) .^ (- (nu + size(x, 2)) / 2), 2);
            y = reshape(y, Nx, Nx);
        end

        function y = studentt_waveform(x, h, mu, sigma, nu)
            % if nargin == 1
            %     A = 600;
            %     sigma = [9, 9];
            %     mu = [10, 10]; % center
            %     nu = 3.2; % degrees of freedom
            % end
            % recover covariance matrix from sigma to Sigma, Sigma is 2*2 matrix
            % Sigma = [sigma(1), conj(sigma(3)); sigma(3), sigma(2)];
            Nd = size(mu, 2);

            if isvector(sigma)
                sigma = sigma(:);

                if length(sigma) == Nd * (Nd + 1) / 2
                    % sigma is tril of Sigma
                    Sigma = tril2full(sigma(:), true, true);
                else
                    % same for the diag
                    Sigma = tril2full(vec2tril(sigma(2:end), -1, Nd), false, true);
                    Sigma = Sigma + diag(sigma(1) .* ones(Nd, 1));
                end

            else
                Sigma = sigma;
            end

            if ndims(x) ~= Nd
                Nx = size(x);
                Nx = Nx(1:end - 1);
                x = reshape(x, [], Nd);
            end

            % calculate Mahalanobis distance

            D = pdist2(x, mu, 'mahalanobis', Sigma);

            % calculate kernel for each data point
            % phi = A .* (1 + D.^2 / nu).^(-(nu + length(f)) / 2);%y = sum(h .* exp(-0.5 * D.^2), 2);
            y = sum(h .* ((1 + D .^ 2) .^ (-nu)), 2);
            y = reshape(y, [Nx(:)', 1]);
        end

        function [y, ycomp] = harmonic_t_flattop(x, h, mu, sigma, nu, d, b, ks, kh, idx_sigma, no_zero_peak, df, taperinfo)

            mu_in = mu;
            sigma_in = sigma;
            h_in = h;
            idx_sigma_in = idx_sigma;

            if nargin < 8 || isempty(idx_sigma)
                idx_sigma = [];
            end

            if nargin < 9 || isempty(no_zero_peak)
                no_zero_peak = false;
            end

            Nk = size(mu, 1); % number of different centers of kernels

            if Nk < length(h)

                for ik = Nk:-1:1

                    if mu(ik) ~= 0
                        % mu_fundmental=mu(ik);
                        mus{ik} = get_harmonic_mu_sigma(mu(ik), 0, ks, kh, no_zero_peak);
                    else
                        mus{ik} = mu(ik);
                    end

                end

                mu = cell2mat(mus)';
            end

            Nd = size(sigma, 2); % dimension of the input data

            Nslice = size(h, Nd + 1); % number of slices of the kernel, used for real and imaginary part of the kernel

            if ndims(x) ~= Nd && size(x, 2) ~= 1
                Nx = size(x);
                Nx = Nx(1:end - 1);
                x = reshape(x, [], Nd);
            else
                Nx = [];
            end

            N = size(x, 1);

            % remove zero height kernel
            if Nd == 2
                h = permute(h, [2, 1, 3]);
                [mux, muy] = meshgrid(mu, mu);
                mask = logical(sum(abs(h) ~= 0, 3));

                mux = mux(mask);
                muy = muy(mask);

                % [found]=ismember(mux+muy,mu);% idx_xy(~found)=[];
                found = true(size(mux, 1), 1);
                mu = [mux, muy];
                mu(~found, :) = [];

                if ~isempty(idx_sigma)
                    idx_sigma = idx_sigma(mask);
                    idx_sigma(~found) = [];
                end

                h = h(logical(mask .* ones([ones(1, Nd), Nslice])));
                h = reshape(h, [], Nslice);
                h(~found, :) = [];
            end

            % remove kernel outside the region

            idx_outside = sum(mu, 2) > max(sum(x, 2));
            % if exist('mu_fundmental','var')
            %     mu_id=round(mu./mu_fundmental);
            %     idx_outside=idx_outside|sum(mu_id, 2)>kh;
            % end
            if size(sigma, 1) == size(mu, 1)
                sigma(idx_outside, :) = [];
            end

            if size(nu, 1) == size(mu, 1)
                nu(idx_outside, :) = [];
            end

            if size(d, 1) == size(mu, 1)
                d(idx_outside, :) = [];
            end

            if size(b, 1) == size(mu, 1)
                b(idx_outside, :) = [];
            end

            mu(idx_outside, :) = [];
            h(idx_outside, :) = [];

            if ~isempty(idx_sigma) && any(idx_outside)
                idx_sigma(idx_outside, :) = [];
            end

            if isempty(h)
                sz = size(x);
                ycomp = zeros([sz(1:(Nd - 1)), Nk]);
                y = zeros([sz(1:(Nd - 1)), 1]);
                return
            end

            % if size(sigma, 2) == Nd * (Nd + 1) / 2
            %     Sigma = tril2full(sigma', true, true);
            % else
            %
            %     Sigma = tril2full(vec2tril((sigma(:, 1) .* sigma(:, 2:end))', -1, Nd), false, true);
            %     Sigma = Sigma + diag2full(sigma(:, 1)' .* ones(Nd, 1));
            % end
            % Sigma=(sigma(:, 1) .* sigma(:, 2:end))';
            Sigma = sigma;
            Sigma(:, 2:end) = sigma(:, 1) .* sigma(:, 2:end);

            if Nd == 2
                Sigma = [Sigma(:, 1), Sigma];
                nu = nu * ones(1, 3);
                d = d * ones(1, 3);
                b = b * ones(1, 3);
            end

            % if any(idx_sigma>Nd)
            %     Sigma=[Sigma;[Sigma(1,1),Sigma(2,2);Sigma(2,2),Sigma(1,1)]];
            %     nu=[[nu,nu];[nu(1),nu(2);nu(2),nu(1)]];
            %     d=[[d,d];[d(1),d(2);d(2),d(1)]];
            %     b=[[b,b];[b(1),b(2);b(2),b(1)]];
            % end

            % if any(idx_sigma>Nd)
            %     Sigma=[Sigma;[Sigma(1,1),Sigma(2,2);Sigma(2,2),Sigma(1,1)]];
            %     nu=[[nu,nu];[nu(1),nu(2);nu(2),nu(1)]];
            %     d=[[d,d];[d(1),d(2);d(2),d(1)]];
            %     b=[[b,b];[b(1),b(2);b(2),b(1)]];
            % end

            % if update_dist
            % D = NaN(N, size(mu, 1));

            h = permute(h, [3, 1, 2]);
            % y = zeros(N, 1, Nslice);
            % h = permute(h, [1, 3, 2]);

            if ~isempty(idx_sigma)
                % Sigma = Sigma(:, :, idx_sigma);
                Sigma = Sigma(idx_sigma, :);
                nu = nu(idx_sigma, :);
                d = d(idx_sigma, :);
                b = b(idx_sigma, :);
            end

            if nargout > 1
                [y, ycomp] = kernel_harmonic_t_flattop(x, h, mu, Sigma, nu, d, b);
                % [y, ycomp] = kernel_harmonic_t_flattop_frf(x, h, mu, Sigma, nu, d, b);
            else
                [y] = kernel_harmonic_t_flattop(x, h, mu, Sigma, nu, d, b);
                % [y] = kernel_harmonic_t_flattop_frf(x, h, mu, Sigma, nu, d, b);
            end

            if size(y, 3) > 1
                y = y(:, :, 1) + 1i .* y(:, :, 2);
            end

            if ~isempty(Nx)
                y = reshape(y, [Nx(:)', 1]);
                y = tril2full(y, false, false);

            end

        end

        function [y, ycomp] = harmonic_t_free(x, h, mu, sigma, nu, d, b, ks, kh, idx_sigma, no_zero_peak, df, taperinfo)

            % mu_in = mu;
            % sigma_in = sigma;
            % h_in = h;
            % idx_sigma_in = idx_sigma;



            Nd = size(sigma, 2); % dimension of the input data

            Nslice = size(h, Nd + 1); % number of slices of the kernel, used for real and imaginary part of the kernel

            if ndims(x) ~= Nd && size(x, 2) ~= 1
                Nx = size(x);
                Nx = Nx(1:end - 1);
                x = reshape(x, [], Nd);
            else
                Nx = [];
            end

            N = size(x, 1);

            % remove zero height kernel
            if Nd == 2
                % h = permute(h, [2, 1, 3]);
                h=reshape(h,[],2);
                idx_nzero=any(h~=0,2);
                h=h(idx_nzero,:);
                mu=mu(idx_nzero,:);
                sigma=sigma(idx_nzero,:);
                nu=nu(idx_nzero,:);
                d=d(idx_nzero,:);
                b=b(idx_nzero,:);
            end





            if isempty(h)
                sz = size(x);
                ycomp = zeros([sz(1:(Nd - 1)), Nk]);
                y = zeros([sz(1:(Nd - 1)), 1]);
                return
            end


            % Sigma = sigma;
            % Sigma(:, 2:end) = sigma(:, 1) .* sigma(:, 2:end);

            if Nd == 2
                % Sigma = [Sigma(:, 1), Sigma];
                % nu = nu * ones(1, 3);
                % d = d * ones(1, 3);
                % b = b * ones(1, 3);
                % sigma=[sigma,b(:,1).*mean(sigma,2)];
                sigma=[sigma,mean(sigma,2)];
                nu=[nu,mean(nu,2)];
                d=[d,mean(d,2)];
                % b=[b,mean(b,2)];
            end



            h = permute(h, [3, 1, 2]);



            if nargout > 1
                [y, ycomp] = kernel_harmonic_t_flattop(x, h, mu, sigma, nu, d, b);
                % [y, ycomp] = kernel_harmonic_t_flattop_frf(x, h, mu, Sigma, nu, d, b);
            else
                [y] = kernel_harmonic_t_flattop(x, h, mu, sigma, nu, d, b);
                % [y] = kernel_harmonic_t_flattop_frf(x, h, mu, Sigma, nu, d, b);
            end

            if size(y, 3) > 1
                y = y(:, :, 1) + 1i .* y(:, :, 2);
            end

            if ~isempty(Nx)
                y = reshape(y, [Nx(:)', 1]);
                y = tril2full(y, false, false);

            end

        end
    end

end
