% Author: Ying Wang, Min Li
% Create Time: 2025
% Copyright (c): 2020-2025 Ying Wang, yingwangrigel@gmail.com,
%                Min Li, minli.231314@gmail.com
% Joint China-Cuba LAB, UESTC, Hangzhou Dianzi University
% License: GNU General Public License v3.0 (see LICENSE file)

% extract the parameter back from vectorized parameter
function [ks, kh, para_xi, para_alpha, para_bs] = vecpara2para(x, order, para_xi, para_alpha, para_bs, model,peak_relation)
if strcmpi(peak_relation, 'harmonic')
    x = x(:);

    if nargin > 5 && ~isempty(model)
        return_stat = true;
    else
        return_stat = false;
    end

    if strcmpi(order, '1+2') || strcmpi(order, '1')
        ks = round(x(1));
        kh = round(x(2));

        para_xi.kernel.vpara = x(2 + 1: ...
            2 + para_xi.kernel.para_info.vecLength);
        para_alpha.kernel.vpara = x(2 + para_xi.kernel.para_info.vecLength + 1: ...
            2 + para_xi.kernel.para_info.vecLength + para_alpha.kernel.para_info.vecLength);

        if ~isempty(para_xi.kernel.fix)
            para_xi.kernel.vpara(~isnan(para_xi.kernel.vfix)) = para_xi.kernel.vfix(~isnan(para_xi.kernel.vfix));
        end

        if ~isempty(para_alpha.kernel.fix)
            para_alpha.kernel.vpara(~isnan(para_alpha.kernel.vfix)) = para_alpha.kernel.vfix(~isnan(para_alpha.kernel.vfix));
        end

        if return_stat
            para_xi.kernel.vt = model.t(2 + 1: ...
                2 + para_xi.kernel.para_info.vecLength);
            para_alpha.kernel.vt = model.t(2 + para_xi.kernel.para_info.vecLength + 1: ...
                2 + para_xi.kernel.para_info.vecLength + para_alpha.kernel.para_info.vecLength);

            % para_xi.kernel.vres = model.res(2 + 1: ...
            %     2 + para_xi.kernel.para_info.vecLength);
            % para_alpha.kernel.vres = model.res(2 + para_xi.kernel.para_info.vecLength + 1: ...
            %     2 + para_xi.kernel.para_info.vecLength + para_alpha.kernel.para_info.vecLength);

            % para_xi.kernel.vjac = model.jac(2 + 1: ...
            %     2 + para_xi.kernel.para_info.vecLength, :);
            % para_alpha.kernel.vjac = model.jac(2 + para_xi.kernel.para_info.vecLength + 1: ...
            %     2 + para_xi.kernel.para_info.vecLength + para_alpha.kernel.para_info.vecLength, :);

        end

    end

    if strcmpi(order, '1+2')
        ks = para_alpha.kernel.info.ks;
        kh = para_alpha.kernel.info.kh;
        para_bs.kernel.vpara = x(2 + para_xi.kernel.para_info.vecLength + para_alpha.kernel.para_info.vecLength + 1: ...
            2 + para_xi.kernel.para_info.vecLength + para_alpha.kernel.para_info.vecLength + para_bs.kernel.para_info.vecLength);

        if ~isempty(para_bs.kernel.fix)
            para_bs.kernel.vpara(~isnan(para_bs.kernel.vfix)) = para_bs.kernel.vfix(~isnan(para_bs.kernel.vfix));
        end

        if return_stat
            para_bs.kernel.vt = model.t(2 + para_xi.kernel.para_info.vecLength + para_alpha.kernel.para_info.vecLength + 1: ...
                2 + para_xi.kernel.para_info.vecLength + para_alpha.kernel.para_info.vecLength + para_bs.kernel.para_info.vecLength);
            % para_bs.kernel.vres = model.res(2 + para_xi.kernel.para_info.vecLength + para_alpha.kernel.para_info.vecLength + 1: ...
            %     2 + para_xi.kernel.para_info.vecLength + para_alpha.kernel.para_info.vecLength + para_bs.kernel.para_info.vecLength);
            % para_bs.kernel.vjac = model.jac(2 + para_xi.kernel.para_info.vecLength + para_alpha.kernel.para_info.vecLength + 1: ...
            %     2 + para_xi.kernel.para_info.vecLength + para_alpha.kernel.para_info.vecLength + para_bs.kernel.para_info.vecLength, :);
        end

    end

    if strcmpi(order, '2')
        ks = para_alpha.kernel.info.ks;
        kh = para_alpha.kernel.info.kh;
        para_bs.kernel.vpara = x(1: ...
            para_bs.kernel.para_info.vecLength);

        if ~isempty(para_bs.kernel.fix)
            para_bs.kernel.vpara(~isnan(para_bs.kernel.vfix)) = para_bs.kernel.vfix(~isnan(para_bs.kernel.vfix));
        end

        if return_stat
            para_bs.kernel.vt = model.t(1: ...
                para_bs.kernel.para_info.vecLength);
            % para_bs.kernel.vres = model.res(1: ...
            %     para_bs.kernel.para_info.vecLength);
            % para_bs.kernel.vjac = model.jac(1: ...
            %     para_bs.kernel.para_info.vecLength, :);
        end

    end

    if strcmpi(order, '1+2') || strcmpi(order, '2')

        if ~para_bs.use_bssigma
            para_bs.kernel.para.sigma(1, 1) = [para_xi.kernel.para.sigma];
            para_bs.kernel.para.sigma(2, 1) = [para_alpha.kernel.para.sigma];
            para_bs.kernel.lb.sigma = para_bs.kernel.para.sigma;
            para_bs.kernel.ub.sigma = para_bs.kernel.para.sigma;
            % else
            %     para_bs.kernel.lb.sigma(1, 1) = [para_xi.kernel.lb.sigma]; %* 2
            %     para_bs.kernel.ub.sigma(1, 1) = [para_xi.kernel.ub.sigma];
            %     para_bs.kernel.lb.sigma(2, 1) = [para_alpha.kernel.lb.sigma];
            %     para_bs.kernel.ub.sigma(2, 1) = [para_alpha.kernel.ub.sigma];
        end

        if ~para_bs.use_bsnu
            para_bs.kernel.para.nu(1) = [para_xi.kernel.para.nu]; %* 2
            para_bs.kernel.para.nu(2) = [para_alpha.kernel.para.nu]; %* 2
            para_bs.kernel.lb.nu = para_bs.kernel.para.nu;
            para_bs.kernel.ub.nu = para_bs.kernel.para.nu;
            % else
            %     para_bs.kernel.lb.nu(1) = para_xi.kernel.lb.nu(1);
            %     para_bs.kernel.ub.nu(1) = para_xi.kernel.ub.nu(1);
            %     para_bs.kernel.lb.nu(2) = para_alpha.kernel.lb.nu(1);
            %     para_bs.kernel.ub.nu(2) = para_alpha.kernel.ub.nu(1);
        end

        if ~para_xi.use_xid
            para_xi.kernel.para.d = para_alpha.kernel.para.d;
            para_xi.kernel.lb.d = para_alpha.kernel.para.d;
            para_xi.kernel.ub.d = para_alpha.kernel.para.d;
            % else
            %     para_xi.kernel.lb.d(1) = para_alpha.kernel.lb.d;
            %     para_xi.kernel.ub.d(1) = para_alpha.kernel.ub.d;
        end

        if ~para_bs.use_bsd
            para_bs.kernel.para.d(1) = [para_xi.kernel.para.d];
            para_bs.kernel.para.d(2) = [para_alpha.kernel.para.d];
            para_bs.kernel.lb.d = para_bs.kernel.para.d;
            para_bs.kernel.ub.d = para_bs.kernel.para.d;
            % else
            %     para_bs.kernel.lb.d(1) = para_xi.kernel.lb.d;
            %     para_bs.kernel.ub.d(1) = para_xi.kernel.ub.d;
            %     para_bs.kernel.lb.d(2) = para_alpha.kernel.lb.d;
            %     para_bs.kernel.ub.d(2) = para_alpha.kernel.ub.d;
        end

        if ~para_bs.use_bsb
            para_bs.kernel.para.b(1) = [para_xi.kernel.para.b];
            para_bs.kernel.para.b(2) = [para_alpha.kernel.para.b];
            para_bs.kernel.lb.b = para_bs.kernel.para.b;
            para_bs.kernel.ub.b = para_bs.kernel.para.b;
            % else
            %     para_bs.kernel.lb.b(1) = para_xi.kernel.lb.b;
            %     para_bs.kernel.ub.b(1) = para_xi.kernel.ub.b;
            %     para_bs.kernel.lb.b(2) = para_alpha.kernel.lb.b;
            %     para_bs.kernel.ub.b(2) = para_alpha.kernel.ub.b;
        end

        para_bs.kernel.lb.mu(2) = para_alpha.kernel.para.mu;
        para_bs.kernel.ub.mu(2) = para_alpha.kernel.para.mu;
        para_bs.kernel.para.mu(2) = para_alpha.kernel.para.mu;
        % para_bs.kernel.para.h=tril2full(para_bs.kernel.para.h,false,false);
        % para_bs.kernel.para.h=permute(para_bs.kernel.para.h,[2,1,3]);

        % para_bs.kernel.para.h(:, :, 1) = tril2full(para_bs.kernel.para.h(:, :, 1), false, false);
        % para_bs.kernel.para.h(:, :, 2) = tril2full(para_bs.kernel.para.h(:, :, 2), false, false);
    end
elseif strcmpi(peak_relation, 'free')
    x = x(:);

    if nargin > 5 && ~isempty(model)
        return_stat = true;
    else
        return_stat = false;
    end

    if strcmpi(order, '1+2') || strcmpi(order, '1')
        ks = round(x(1));
        kh = round(x(2));

        para_xi.kernel.vpara = x(2 + 1: ...
            2 + para_xi.kernel.para_info.vecLength);
        para_alpha.kernel.vpara = x(2 + para_xi.kernel.para_info.vecLength + 1: ...
            2 + para_xi.kernel.para_info.vecLength + para_alpha.kernel.para_info.vecLength);

        if ~isempty(para_xi.kernel.fix)
            para_xi.kernel.vpara(~isnan(para_xi.kernel.vfix)) = para_xi.kernel.vfix(~isnan(para_xi.kernel.vfix));
        end

        if ~isempty(para_alpha.kernel.fix)
            para_alpha.kernel.vpara(~isnan(para_alpha.kernel.vfix)) = para_alpha.kernel.vfix(~isnan(para_alpha.kernel.vfix));
        end

        if return_stat
            para_xi.kernel.vt = model.t(2 + 1: ...
                2 + para_xi.kernel.para_info.vecLength);
            para_alpha.kernel.vt = model.t(2 + para_xi.kernel.para_info.vecLength + 1: ...
                2 + para_xi.kernel.para_info.vecLength + para_alpha.kernel.para_info.vecLength);

            % para_xi.kernel.vres = model.res(2 + 1: ...
            %     2 + para_xi.kernel.para_info.vecLength);
            % para_alpha.kernel.vres = model.res(2 + para_xi.kernel.para_info.vecLength + 1: ...
            %     2 + para_xi.kernel.para_info.vecLength + para_alpha.kernel.para_info.vecLength);

            % para_xi.kernel.vjac = model.jac(2 + 1: ...
            %     2 + para_xi.kernel.para_info.vecLength, :);
            % para_alpha.kernel.vjac = model.jac(2 + para_xi.kernel.para_info.vecLength + 1: ...
            %     2 + para_xi.kernel.para_info.vecLength + para_alpha.kernel.para_info.vecLength, :);

        end

    end

    if strcmpi(order, '1+2')
        ks = para_alpha.kernel.info.ks;
        kh = para_alpha.kernel.info.kh;
        para_bs.kernel.vpara = x(2 + para_xi.kernel.para_info.vecLength + para_alpha.kernel.para_info.vecLength + 1: ...
            2 + para_xi.kernel.para_info.vecLength + para_alpha.kernel.para_info.vecLength + para_bs.kernel.para_info.vecLength);

        if ~isempty(para_bs.kernel.fix)
            para_bs.kernel.vpara(~isnan(para_bs.kernel.vfix)) = para_bs.kernel.vfix(~isnan(para_bs.kernel.vfix));
        end

        if return_stat
            para_bs.kernel.vt = model.t(2 + para_xi.kernel.para_info.vecLength + para_alpha.kernel.para_info.vecLength + 1: ...
                2 + para_xi.kernel.para_info.vecLength + para_alpha.kernel.para_info.vecLength + para_bs.kernel.para_info.vecLength);
            % para_bs.kernel.vres = model.res(2 + para_xi.kernel.para_info.vecLength + para_alpha.kernel.para_info.vecLength + 1: ...
            %     2 + para_xi.kernel.para_info.vecLength + para_alpha.kernel.para_info.vecLength + para_bs.kernel.para_info.vecLength);
            % para_bs.kernel.vjac = model.jac(2 + para_xi.kernel.para_info.vecLength + para_alpha.kernel.para_info.vecLength + 1: ...
            %     2 + para_xi.kernel.para_info.vecLength + para_alpha.kernel.para_info.vecLength + para_bs.kernel.para_info.vecLength, :);
        end

    end

    if strcmpi(order, '2')
        ks = para_alpha.kernel.info.ks;
        kh = para_alpha.kernel.info.kh;
        para_bs.kernel.vpara = x(1: ...
            para_bs.kernel.para_info.vecLength);

        if ~isempty(para_bs.kernel.fix)
            para_bs.kernel.vpara(~isnan(para_bs.kernel.vfix)) = para_bs.kernel.vfix(~isnan(para_bs.kernel.vfix));
        end

        if return_stat
            para_bs.kernel.vt = model.t(1: ...
                para_bs.kernel.para_info.vecLength);
            % para_bs.kernel.vres = model.res(1: ...
            %     para_bs.kernel.para_info.vecLength);
            % para_bs.kernel.vjac = model.jac(1: ...
            %     para_bs.kernel.para_info.vecLength, :);
        end

    end

    if strcmpi(order, '1+2') || strcmpi(order, '2')

        if ~para_bs.use_bssigma
            % para_bs.kernel.para.sigma(1, 1) = [para_xi.kernel.para.sigma];
            % para_bs.kernel.para.sigma(2, 1) = [para_alpha.kernel.para.sigma];

            sigma1=[para_xi.kernel.para.sigma; para_alpha.kernel.para.sigma];
            para_bs.kernel.para.sigma=get_ndgrid([sigma1,sigma1], 'list');

            para_bs.kernel.lb.sigma = para_bs.kernel.para.sigma;
            para_bs.kernel.ub.sigma = para_bs.kernel.para.sigma;
            % else
            %     para_bs.kernel.lb.sigma(1, 1) = [para_xi.kernel.lb.sigma]; %* 2
            %     para_bs.kernel.ub.sigma(1, 1) = [para_xi.kernel.ub.sigma];
            %     para_bs.kernel.lb.sigma(2, 1) = [para_alpha.kernel.lb.sigma];
            %     para_bs.kernel.ub.sigma(2, 1) = [para_alpha.kernel.ub.sigma];
        else
            idx_xi=[true(size(para_xi.kernel.para.sigma)); false(size(para_alpha.kernel.para.sigma))];
            idx_xi=get_ndgrid([idx_xi,idx_xi], 'list');
            sigma1=[para_xi.kernel.para.sigma; para_alpha.kernel.para.sigma];
            sigma=get_ndgrid([sigma1,sigma1], 'list');
            para_bs.kernel.para.sigma(idx_xi)=sigma(idx_xi);

            para_bs.kernel.lb.sigma = para_bs.kernel.para.sigma;
            para_bs.kernel.ub.sigma = para_bs.kernel.para.sigma;
        end

        if ~para_bs.use_bsnu
            % para_bs.kernel.para.nu(1) = [para_xi.kernel.para.nu]; %* 2
            % para_bs.kernel.para.nu(2) = [para_alpha.kernel.para.nu]; %* 2
            nu1=[para_xi.kernel.para.nu; para_alpha.kernel.para.nu];
            para_bs.kernel.para.nu=get_ndgrid([nu1,nu1], 'list');

            para_bs.kernel.lb.nu = para_bs.kernel.para.nu;
            para_bs.kernel.ub.nu = para_bs.kernel.para.nu;
            % else
            %     para_bs.kernel.lb.nu(1) = para_xi.kernel.lb.nu(1);
            %     para_bs.kernel.ub.nu(1) = para_xi.kernel.ub.nu(1);
            %     para_bs.kernel.lb.nu(2) = para_alpha.kernel.lb.nu(1);
            %     para_bs.kernel.ub.nu(2) = para_alpha.kernel.ub.nu(1);
        else
            idx_xi=[true(size(para_xi.kernel.para.nu)); false(size(para_alpha.kernel.para.nu))];
            idx_xi=get_ndgrid([idx_xi,idx_xi], 'list');
            nu1=[para_xi.kernel.para.nu; para_alpha.kernel.para.nu];
            nu=get_ndgrid([nu1,nu1], 'list');
            para_bs.kernel.para.nu(idx_xi)=nu(idx_xi);
        end

        if ~para_xi.use_xid
            para_xi.kernel.para.d = para_alpha.kernel.para.d;
            para_xi.kernel.lb.d = para_alpha.kernel.para.d;
            para_xi.kernel.ub.d = para_alpha.kernel.para.d;
            % else
            %     para_xi.kernel.lb.d(1) = para_alpha.kernel.lb.d;
            %     para_xi.kernel.ub.d(1) = para_alpha.kernel.ub.d;
        end

        if ~para_bs.use_bsd
            % para_bs.kernel.para.d(1) = [para_xi.kernel.para.d];
            % para_bs.kernel.para.d(2) = [para_alpha.kernel.para.d];

            d1=[para_xi.kernel.para.d; para_alpha.kernel.para.d];
            para_bs.kernel.para.d=get_ndgrid([d1,d1], 'list');

            para_bs.kernel.lb.d = para_bs.kernel.para.d;
            para_bs.kernel.ub.d = para_bs.kernel.para.d;
            % else
            %     para_bs.kernel.lb.d(1) = para_xi.kernel.lb.d;
            %     para_bs.kernel.ub.d(1) = para_xi.kernel.ub.d;
            %     para_bs.kernel.lb.d(2) = para_alpha.kernel.lb.d;
            %     para_bs.kernel.ub.d(2) = para_alpha.kernel.ub.d;
        else
            idx_xi=[true(size(para_xi.kernel.para.d)); false(size(para_alpha.kernel.para.d))];
            idx_xi=get_ndgrid([idx_xi,idx_xi], 'list');
            d1=[para_xi.kernel.para.d; para_alpha.kernel.para.d];
            d=get_ndgrid([d1,d1], 'list');
            para_bs.kernel.para.d(idx_xi)=d(idx_xi);
        end

        if ~para_bs.use_bsb
            % para_bs.kernel.para.b(1) = [para_xi.kernel.para.b];
            % para_bs.kernel.para.b(2) = [para_alpha.kernel.para.b];
            b1=[para_xi.kernel.para.b; para_alpha.kernel.para.b];
            para_bs.kernel.para.b=get_ndgrid([b1,b1], 'list');

            para_bs.kernel.lb.b = para_bs.kernel.para.b;
            para_bs.kernel.ub.b = para_bs.kernel.para.b;
            % else
            %     para_bs.kernel.lb.b(1) = para_xi.kernel.lb.b;
            %     para_bs.kernel.ub.b(1) = para_xi.kernel.ub.b;
            %     para_bs.kernel.lb.b(2) = para_alpha.kernel.lb.b;
            %     para_bs.kernel.ub.b(2) = para_alpha.kernel.ub.b;
        end
        mu1=[para_xi.kernel.para.mu; para_alpha.kernel.para.mu];
        para_bs.kernel.para.mu = get_ndgrid([mu1,mu1], 'list');
        para_bs.kernel.lb.mu = para_bs.kernel.para.mu;
        para_bs.kernel.ub.mu = para_bs.kernel.para.mu;

        % para_bs.kernel.lb.mu(2) = para_alpha.kernel.para.mu;
        % para_bs.kernel.ub.mu(2) = para_alpha.kernel.para.mu;
        % para_bs.kernel.para.mu(2) = para_alpha.kernel.para.mu;
        % para_bs.kernel.para.h=tril2full(para_bs.kernel.para.h,false,false);
        % para_bs.kernel.para.h=permute(para_bs.kernel.para.h,[2,1,3]);

        % para_bs.kernel.para.h(:, :, 1) = tril2full(para_bs.kernel.para.h(:, :, 1), false, false);
        % para_bs.kernel.para.h(:, :, 2) = tril2full(para_bs.kernel.para.h(:, :, 2), false, false);
    end
end
end

% ksold=ks;
% khold=kh;

% if strcmpi(order, '1+2') || strcmpi(order, '1')

%     ks = round(x(1));
%     ks2max = ks-ksmax;
%     kh = round(x(2));% kh = ceil(x(1));
%     kh2max = kh-khmax;
%     % k=kh+ks;
%     kmax=khmax+ksmax;

%     if kh2max > 0
%         % warning('kh is reching the maximum')
%         % kh2max=0;
%         kh=khold;
%     end
%     if ks2max > 0
%         % warning('ks is reching the maximum')
%         % ks2max=0;
%         ks=ksold;
%     end

%     if strcmpi(para_xi.type, 'tstudent')

%         para_xi.kernel.para.h = x(3);
%         para_xi.kernel.para.sigma = x(4);
%         para_xi.kernel.para.nu = x(5);
%         para_xi.kernel.para.mu = x(6);
%         para_alpha.kernel.para.mu = x(7);
%         para_alpha.kernel.para.sigma = x(8);
%         para_alpha.kernel.alpha = x(9:9 + kmax - 1).';

%         % if kh2max < 0
%         %     para_alpha.kernel.alpha = para_alpha.kernel.alpha(1:ksold+kh);
%         % end
%         % if ks2max < 0
%         %     para_alpha.kernel.alpha(1:-ks2max)=[];
%         % end

%         L1 = 9 + kmax;
%     elseif strcmpi(para_xi.type, 'tstudent+peak')

%         para_xi.kernel.para.h = x(3);
%         para_xi.kernel.para.sigma = x(4);
%         para_xi.kernel.para.nu = x(5);
%         para_xi.kernel.para.mu = x(6);
%         para_xi.kernel.para.mu = x(7);
%         para_xi.kernel.para.sigma = x(8);
%         para_xi.alpha = x(9);
%         para_alpha.kernel.para.mu = x(10);
%         para_alpha.kernel.para.sigma = x(11);
%         para_alpha.kernel.alpha = x(12:12 + kmax - 1).';

%         % if kh2max < 0
%         %     para_alpha.kernel.alpha = para_alpha.kernel.alpha(1:ksold+kh);
%         % end
%         % if ks2max < 0
%         %     para_alpha.kernel.alpha(1:-ks2max)=[];
%         % end

%         L1 = 12 + kmax;
%     elseif strcmpi(para_xi.type, 'gaussian')
%         para_xi.kernel.vpara=x(3:5);
%         para_alpha.kernel.para.mu = x(6);
%         para_alpha.kernel.para.sigma = x(7);
%         para_alpha.kernel.alpha = x(8:8 + kmax - 1).';
%         L1 = 8 + kmax;

%     end
%     idx=get_alpha_idx(0,ksmax,khmax,ks,kh);
%     para_alpha.kernel.alpha = para_alpha.kernel.alpha(idx);
% else
%     L1 = 1;
%     % para_xi = [];
%     % para_alpha = [];
% end

% if strcmpi(order, '1+2') || strcmpi(order, '2')
%     % this step doesn't change kh and ks
%     % kh2max = kh-khmax;
%     % ks2max = ks-ksmax;
%     % k=kh+ks;
%     kmax=khmax+ksmax;
%     % N1=k+kxi;
%     N1MAX=kmax+kxi;

%     % N2 = N1 ^ 2;
%     N2MAX = N1MAX ^ 2;

%     para_bs.kernel.para.h(:,:,1) = x(L1:L1 + N2MAX - 1);
%     para_bs.kernel.para.h(:,:,2) = x(L1 + N2MAX:L1 + N2MAX + N2MAX - 1);
%     para_bs.kernel.para.h(:,:,1) = reshape(para_bs.kernel.para.h(:,:,1), [N1MAX, N1MAX]);
%     para_bs.kernel.para.h(:,:,2) = reshape(para_bs.kernel.para.h(:,:,2), [N1MAX, N1MAX]);
%     para_bs.betaxy = reshape(x(L1 + N2MAX + N2MAX:L1 + N2MAX + N2MAX +N1MAX - 1), 1, []);
%     para_bs.gammaxy = reshape(x(L1 + N2MAX + N2MAX + N1MAX:L1 + N2MAX + N2MAX + 2 * N1MAX - 1), 1, []);

%     if para_bs.use_bssigma
%         if strcmpi(para_xi.type, 'tstudent+peak')
%             para_bs.xi_B = x(end - 2);
%             para_bs.xi_sigma = x(end - 1);
%             para_bs.alpha_sigma = x(end);
%             if L1 + N2MAX + N2MAX +2 * N1MAX - 1+3~=length(x)
%                 error('para length not fit')
%             end
%         elseif strcmpi(para_xi.type, 'tstudent')
%             para_bs.xi_B = x(end - 1);
%             para_bs.alpha_sigma = x(end);
%             if L1 + N2MAX + N2MAX +2 * N1MAX - 1+2~=length(x)
%                 error('para length not fit')
%             end
%         elseif strcmpi(para_xi.type, 'gaussian')
%             para_bs.xi_sigma = x(end - 1);
%             para_bs.alpha_sigma = x(end);
%             if L1 + N2MAX + N2MAX +2 * N1MAX - 1+2~=length(x)
%                 error('para length not fit')
%             end
%         else
%             error('no such method')
%         end
%     end
%     idx=get_alpha_idx(kxi,ksmax,khmax,ks,kh);
%     para_bs.kernel.para.h(:,:,1) = para_bs.kernel.para.h(:,:,1)(idx,idx);
%     para_bs.kernel.para.h(:,:,2) = para_bs.kernel.para.h(:,:,2)(idx,idx);
%     para_bs.betaxy = para_bs.betaxy(idx);
%     para_bs.gammaxy = para_bs.gammaxy(idx);
% else
%     para_bs.kernel.para.h(:,:,1) = [];
%     para_bs.kernel.para.h(:,:,2) = [];
%     para_bs.betaxy = [];
%     para_bs.gammaxy = [];
% end

% para_bs.kernel.para.h(:,:,1) = tril2full(para_bs.kernel.para.h(:,:,1));
% para_bs.kernel.para.h(:,:,2) = tril2full(para_bs.kernel.para.h(:,:,2));

% end

%{
% extract the parameter back from vectorized parameter
function [kh, para_xi, para_alpha, para_bs] = vecpara2para(x, kh, order, hos_components, para_xi, para_alpha, para_bs)
x = x(:);

% warning('kh is all from the initial, but not the last loop')

if strcmpi(order, '1+2') || strcmpi(order, '1')

    knew  = round(x(1));
    kdiff=knew-kh;

    if strcmpi(para_xi.type, 'tstudent')

        para_xi.kernel.para.h = x(2);
        para_xi.kernel.para.sigma = x(3);
        para_xi.kernel.para.nu = x(4);
        para_xi.kernel.para.mu = x(5);
        para_alpha.kernel.para.mu = x(6);
        para_alpha.kernel.para.sigma = x(7);
        if kdiff>=0
            para_alpha.kernel.alpha = padarray(x(8:8 + kh - 1).',[0,kdiff],"replicate",'post');
        else
            para_alpha.kernel.alpha = x(8:8 + kh - 1).';
            para_alpha.kernel.alpha = para_alpha.kernel.alpha(1:knew);
        end
        L1 = 8 + kh;
    elseif strcmpi(para_xi.type, 'tstudent+peak')

        para_xi.kernel.para.h = x(2);
        para_xi.kernel.para.sigma = x(3);
        para_xi.kernel.para.nu = x(4);
        para_xi.kernel.para.mu = x(5);
        para_xi.kernel.para.mu = x(6);
        para_xi.kernel.para.sigma = x(7);
        para_xi.alpha = x(8);
        para_alpha.kernel.para.mu = x(9);
        para_alpha.kernel.para.sigma = x(10);

        if kdiff>0
            para_alpha.kernel.alpha = padarray(x(11:11 + kh - 1).',[0,kdiff],"replicate",'post');
        elseif kdiff<0
            para_alpha.kernel.alpha = x(11:11 + kh - 1).';
            para_alpha.kernel.alpha = para_alpha.kernel.alpha(1:knew);
        end
        L1 = 11 + kh;
    end

else
    L1 = 1;
    para_xi = [];
    para_alpha = [];
end

if strcmpi(order, '1+2') || strcmpi(order, '2')

    if hos_components

        if strcmpi(para_bs.type, 'bigaussian+xi+delta')
            N1=kh+2;
        elseif strcmpi(para_bs.type, 'bigaussian+xi')||strcmpi(para_bs.type, 'trigaussian+xi')
            N1 = kh + 1;
        end

    else
        N1 = kh;
    end
    N2 = N1 ^ 2;
    para_bs.kernel.para.h(:,:,1) = x(L1:L1 + N2 - 1);
    para_bs.kernel.para.h(:,:,2) = x(L1 + N2:L1 + N2 + N2 - 1);
    para_bs.kernel.para.h(:,:,1) = reshape(para_bs.kernel.para.h(:,:,1), [N1, N1]);
    para_bs.kernel.para.h(:,:,2) = reshape(para_bs.kernel.para.h(:,:,2), [N1, N1]);
    para_bs.betaxy = reshape(x(L1 + N2 + N2:L1 + N2 + N2 +N1 - 1), 1, []);
    para_bs.gammaxy = reshape(x(L1 + N2 + N2 + N1:L1 + N2 + N2 +2 * N1 - 1), 1, []);
    para_bs.xi_B = x(end-2);
    para_bs.xi_sigma= x(end-1);
    para_bs.alpha_sigma = x(end);

    if exist('kdiff','var')
        if kdiff>0
            para_bs.kernel.para.h(:,:,1) = padarray(para_bs.kernel.para.h(:,:,1),[kdiff,kdiff],"replicate",'post');
            para_bs.kernel.para.h(:,:,2) = padarray(para_bs.kernel.para.h(:,:,2),[kdiff,kdiff],"replicate",'post');
            para_bs.betaxy = padarray(para_bs.betaxy,[0,kdiff],"replicate",'post');
            para_bs.gammaxy = padarray(para_bs.gammaxy,[0,kdiff],"replicate",'post');
        elseif kdiff<0
            para_bs.kernel.para.h(:,:,1) = para_bs.kernel.para.h(:,:,1)(1:end+kdiff,1:end+kdiff);
            para_bs.kernel.para.h(:,:,2) = para_bs.kernel.para.h(:,:,2)(1:end+kdiff,1:end+kdiff);
            para_bs.betaxy = para_bs.betaxy(1,1:end+kdiff);
            para_bs.gammaxy = para_bs.gammaxy(1,1:end+kdiff);

        end
    end
else
    para_bs.kernel.para.h(:,:,1) = [];
    para_bs.kernel.para.h(:,:,2) = [];
    para_bs.betaxy=[];
    para_bs.gammaxy=[];
end

para_bs.kernel.para.h(:,:,1)=tril2full(para_bs.kernel.para.h(:,:,1));
para_bs.kernel.para.h(:,:,2)=tril2full(para_bs.kernel.para.h(:,:,2));

if exist('knew','var')
    kh=knew;
end

end

%}
