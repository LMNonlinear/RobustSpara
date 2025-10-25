function loss = ObjFunBiXiAlpha(self, para, record)

if self.iiter == 0
    self.temp.para = NaN(length(para), 2500);
else
    self.temp.para(:, int64(self.iiter + 1)) = para;
end

% try

if nargin < 3
    record = [];
end

% record: avoid redundant calculation
% get para to fit from vectorized para
para = recover_fixed_para(para, [], [], self.para_fit.parainfo);
[ks, kh, para_xi, para_alpha, para_bs] = vecpara2para(para, self.order, self.para_xi, self.para_alpha, self.para_bs,[],self.peak_relation);

%% spectrum
if strcmpi(self.order, '1')
    shat = pred_s(self.f, ks, kh, para_xi, para_alpha);
    % elseif strcmpi(self.order, '2')
    % shat_xi = pred_s_xi(self.f, para_xi);
elseif strcmpi(self.order, '1+2')
    shat = pred_s(self.f, ks, kh, para_xi, para_alpha);
    % comp = pred_comp(self.f, ks, kh, para_xi, para_alpha, para_bs, self.hos_components);
    % shat_xi = pred_s_xi(self.f, para_xi);
end

%% bispectrum
if strcmpi(self.order, '1+2')
    % bs = self.bs;
    bshat = pred_bs(self.f, para_bs, self.fxfy);

    switch self.loss_type
        case {'bicoherence'}

            bshatdenom = sqrt(get_sxsysxy(self.f, shat, [], self.fxy)); %sxsysxyhat
            bshatdenom(bshatdenom == 0) = NaN;

            bc = record.bc;
            bshatvar = bshatdenom;

        otherwise
            % bsvar = self.bsvar;
    end

    if strcmpi(self.region_bs, "mu+diag")

        if self.para_alpha.isubharmonic %&& (strcmpi(self.hos_components,'xi+alpha') || strcmpi(self.hos_components, 'alpha'))
            idx_mu = find_closest(self.f, self.para_alpha.kernel.para.mu * 2);
        else
            idx_mu = find_closest(self.f, self.para_alpha.kernel.para.mu);
        end

        % bs = [bs(idx_mu, :)'; record.bs_diag];
        bshat = [bshat(idx_mu, :)'; bshat(record.idx_diag)];
        % bsvar = [bsvar(idx_mu, :)'; record.bsvar_diag];
        bshatvar = [bshatvar(idx_mu, :)'; record.bshatvar_diag];
    end

elseif strcmpi(self.order, '2')
    bshat = pred_bs(self.f, para_bs, self.fxfy);

    if strcmpi(self.region_bs, "boundary")
        % bs = record.bs_boundary;
        bshat = bshat(record.idx_mu, :)';
        % bsvar = record.bsvar_boundary;
        bshatvar = record.bshatvar_boundary;
    elseif strcmpi(self.region_bs, "diag")
        % bs = record.bs_diag;
        bshat = bshat(record.idx_diag);
        % bsvar = record.bsvar_diag;
        bshatvar = record.bshatvar_diag;
    elseif strcmpi(self.region_bs, "offdiag") || strcmpi(self.region_bs, "all") || isempty(self.region_bs)

        % bs = self.bs;
        % bshat=bshat;
        bc = record.bc;
        % bsvar = record.bsvar;
        bshatvar = record.bshatvar;
    else
        error('not support')
    end

end

%% lost function
% if strcmpi(self.order, '1+2')

%     switch self.loss_type
%         case {'Hetero', 'Heteroscedasticity'}
%             bshatvar = [];
%         otherwise
%             % reweighted, use last iter
%             % if self.iiter>1
%             %     bshatvar_pre=self.para_bs.bshatvar;
%             %     self.para_bs.bshatdenom_pre=self.para_bs.bshatdenom;
%             %     self.para_bs.bshatdenom=bshatdenom;
%             %     self.para_bs.bshatvar=bshatvar;
%             %     bshatvar=bshatvar_pre;
%             % else
%             self.para_bs.bshatdenom = bshatdenom;
%             self.para_bs.bshatvar = bshatvar;
%             % end
%     end

% end

if strcmpi(self.order, '1')
    %%{
    loss = calc_loss_s(self, self.loss_type, self.s, shat, self.svar); %sep_complex
    % loss = [loss; 0.1 * (calc_loss_s(self, self.loss_type, self.trend, shat_xi, self.svar))];
    %%}
    %{
    svar = ones(size(self.svar));
    loss = calc_loss_s(self, self.loss_type, self.s, shat, svar); %sep_complex
    % loss = (calc_loss_s(self, self.loss_type, repmat(self.s, 2, 1), [shat; shat_xi], repmat(self.svar, 2, 1))); %sep_complex
    loss = [loss; 0.5 * (calc_loss_s(self, self.loss_type, self.trend, shat_xi, svar))];
    %}
elseif strcmpi(self.order, '2')
    % loss = calc_loss_bs(self, self.loss_type, bs, bshat, bsvar, bshatvar); %sep_complex
    loss = calc_loss_bs(self, self.loss_type, bc, bshat, [], bshatvar); %sep_complex
elseif strcmpi(self.order, '1+2')
    %%{
    [loss] = calc_loss(self, self.loss_type, self.s, shat, self.svar, bc, bshat, [], bshatvar);
    % loss = [loss; 0.1 * p .* (calc_loss_s(self, self.loss_type, self.trend, shat_xi, self.svar))]; %
    %%}
    %{
    svar = ones(size(self.svar));
    % [loss, p, q] = calc_loss(self, self.loss_type, self.s, shat, self.svar, bs, bshat, bsvar, bshatvar); %sep_complex
    [loss, p, q] = calc_loss(self, self.loss_type, self.s, shat, svar, bc, bshat, [], bshatvar);
    % loss = (calc_loss(self, self.loss_type, repmat(self.s, 2, 1), [shat; shat_xi], repmat(self.svar, 2, 1), bs, bshat, bsvar, bshatvar)); %sep_complex
    loss = [loss; 0.5 * p .* (calc_loss_s(self, self.loss_type, self.trend, shat_xi, svar))]; %
    % loss = [loss;  (calc_loss_s(self, self.loss_type, self.trend, shat_xi, self.svar))];
    %}
end

penalty = calc_penalty(self, para_xi, para_alpha, para_bs);
loss = [loss; penalty];

self.iiter = self.iiter + 1;

% %{
if strcmpi(self.order, '1')
    % ifit = 1;
    itv_show = 50;

    if self.verbose && (mod(self.iiter, itv_show) == 1 || self.iiter == 1)
        s = meanlog10(self.s, 2);
        figure(10);
        plot(self.f,log10([s, shat]))
        % plot(([s, shat]))
        drawnow
        sgtitle(['iter=', num2str(self.iiter)])
        % pause(0.2)
        if self.iiter == 1
            % pause
        end

    end

elseif strcmpi(self.order, '2')
    itv_show = 400;
    % self.iiter
    if self.verbose && (mod(self.iiter, itv_show) == 1 || self.iiter == 1)
        figure(11)

            %{
            y = [log10([(real(bs(:)) .^ 2) ./ (bshatvar(:)), (real(bshat(:)) .^ 2) ./ (bshatvar(:))]); ...
                log10([(imag(bs(:)) .^ 2) ./ (bshatvar(:)), (imag(bshat(:)) .^ 2) ./ (bshatvar(:))])];
            plot([y])

            if max(y, [], 'all') >- 6
                ylim([-6, max(y, [], 'all')])
            end
            pause(0.1)
            drawnow
            %}
        % bchat = real(pred_bs(self.f, para_bs, self.fxfy) ./ (sqrt(self.shat .* self.shat.' .* self.sxyhat)));
        % bchat = real(pred_bs(self.f, para_bs, self.fxfy) ./ (mean(self.bshatdenom, 3)));
        bchat = real(bshat ./ mean(self.bshatdenom, 3));
        subplot(2, 1, 1)
        surfbc(self.fx, self.fy, real(mean(self.bc, 3)));
        clim('auto')
        subplot(2, 1, 2)
        surfbc(self.fx, self.fy, bchat);
        clim('auto')
        % zlim([min(bchat,[],"all"),max(bchat,[],"all")])
        drawnow
        sgtitle(['iter=', num2str(self.iiter)])
        % pause(0.2)
        if self.iiter == 1
            % pause
        end

    end

elseif strcmpi(self.order, '1+2')
    % %{
    % ifit = 3;
    itv_show = 400;
    % self.iiter
    if self.verbose && (mod(self.iiter, itv_show) == 1 || self.iiter == 1)
        figure(12)
        % y = [log10(self.s),log10(shat)];

        % y = [log10(self.s),log10(shat);
        %     log10([(real(bs(:)) .^ 2) ./ (bshatvar(:)), (real(bshat(:)) .^ 2) ./ (bshatvar(:))]); ...
        %     log10([(imag(bs(:)) .^ 2) ./ (bshatvar(:)), (imag(bshat(:)) .^ 2) ./ (bshatvar(:))])];

            %{
            y = [log10(self.s.^2./shat),log10(shat.^2./shat);
                log10([(real(bs(:)) .^ 2) ./ (bshatvar(:)), (real(bshat(:)) .^ 2) ./ (bshatvar(:))]); ...
                log10([(imag(bs(:)) .^ 2) ./ (bshatvar(:)), (imag(bshat(:)) .^ 2) ./ (bshatvar(:))])];

            plot([y])
            if max(y, [], 'all') >- 6
                ylim([-6, max(y, [], 'all')])
            end
            pause(0.1)
            drawnow
            %}
        subplot(3, 1, 1)
        plot(log10([self.s, shat]))

        % bchat = real(pred_bs(self.f, para_bs, self.fxfy) ./ (mean(bsdenom, 3)));

        % subplot(3, 2, 2)

        subplot(3, 2, 3)
        bc = mean(self.bc, 3);
        surfbc(self.fx, self.fy, real(bc));
        clim('auto')

        subplot(3, 2, 4)
        surfbc(self.fx, self.fy, imag(bc));
        clim('auto')

        subplot(3, 2, 5)
        bchat = (bshat ./ mean(self.bshatdenom, 3));
        surfbc(self.fx, self.fy, real(bchat));
        clim('auto')
        % zlim([min(bchat,[],"all"),max(bchat,[],"all")])
        drawnow
        subplot(3, 2, 6)
        surfbc(self.fx, self.fy, imag(bchat));
        clim('auto')
        % zlim([min(bchat,[],"all"),max(bchat,[],"all")])
        drawnow
        sgtitle(['iter=', num2str(self.iiter)])
        % pause(0.2)
        if self.iiter == 1
            % pause
        end

    end

    % %}
end

% %}

% if self.verbose
self.model.loss_history(self.iiter, :) = norm(loss(~isnan(loss)));
% self.model.para_xi.nu(self.iiter, :) = para_xi.kernel.para.nu;

% end

% if strcmpi(self.order, '1+2') && strcmpi(self.loss_type, 'AIC')
%     loss = 2 * loss + 2 * self.Np;
% end

% if ~isdouble(loss)
%     loss = double(loss);
% end
loss = convert_to_double(loss);

% disp(self.iiter)
% catch ME
%     % warning(ME)
%     disp(getReport(ME, 'extended', 'hyperlinks', 'on'))
% end

end
