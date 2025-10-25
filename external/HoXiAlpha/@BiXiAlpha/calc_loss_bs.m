function loss = calc_loss_bs(self, loss_type, bs, bshat, bsvar, bshatvar)

if isempty(loss_type)
    loss_type = self.loos_type;
end

% if size(bs, 3) == 1
%     idx_calc = self.idxused; %~isnan(bs) & ~isnan(bshat); %&bshat~=0
% else
%     idx_calc = self.idxused_segs;
% end
% idx_calc = self.idxused_segs;
idx_calc = self.tfused_segs;

% if isempty(bsvar)
%     warning('bsvar is empty, use homoskedasticity')
%     bsvar = ones(size(bs));
% end

% if nargin > 4 && ~isempty(bshatvar)
%     idx_calc = idx_calc & ~isnan(bshatvar) & ~isnan(bsvar) & bshatvar ~= 0 & bsvar ~= 0;
% end

switch loss_type
    case {'bicoherence'}
        %{
        % bs = bs(idx_calc);
        % bshat = bshat(idx_calc);
        bsstd = sqrt(bsvar) + self.dn_reg; %sqrt
        bshatstd = sqrt(bshatvar) + self.dn_reg; %sqrt
        loss = cat(3, (real(bs) ./ bsstd - real(bshat) ./ bshatstd), ...
            (imag(bs) ./ bsstd - imag(bshat) ./ bshatstd)); %abs
        % loss = sqrt((bchat(idx_calc) - bc(idx_calc)) .^ 2);
        idx_calc = self.tfused_segs_sep_z;
        %}
        bc = bs;
        % TODO: input is not bshatvar, already bshatstd, so bshatdenom = bshatvar
        bshatdenom = bshatvar + self.dn_reg; %sqrt
        % %{
        loss = cat(3, (real(bc) - real(bshat) ./ bshatdenom), ...
            (imag(bc) - imag(bshat) ./ bshatdenom)); %abs
        idx_calc = self.tfused_segs_sep_z;
        % %}
        %{
        % loss = abs(bc-bshat ./ bshatdenom); %abs
        % loss = loss(idx_calc);
        %}
        % loss = cat(3, abs(real(bs) ./ bsdenom - real(bshat) ./ bshatdenom).^2, ...
        %               abs(imag(bs) ./ bsdenom - imag(bshat) ./ bshatdenom).^2); %abs
    case {'bicoherence_complecated'}
        %{
        % bs = bs(idx_calc);
        % bshat = bshat(idx_calc);
        bsstd = sqrt(bsvar) + self.dn_reg; %sqrt
        bshatstd = sqrt(bshatvar) + self.dn_reg; %sqrt
        loss = cat(3, (real(bs) ./ bsstd - real(bshat) ./ bshatstd), ...
            (imag(bs) ./ bsstd - imag(bshat) ./ bshatstd)); %abs
        % loss = sqrt((bchat(idx_calc) - bc(idx_calc)) .^ 2);
        idx_calc = self.tfused_segs_sep_z;
        %}

        bsdenom = bsvar + self.dn_reg; %sqrt
        bshatdenom = bshatvar + self.dn_reg; %sqrt
        loss = cat(3, (real(bs) ./ bsdenom - real(bshat) ./ bshatdenom), ...
            (imag(bs) ./ bsdenom - imag(bshat) ./ bshatdenom)); %abs
        % loss = cat(3, abs(real(bs) ./ bsdenom - real(bshat) ./ bshatdenom).^2, ...
        %               abs(imag(bs) ./ bsdenom - imag(bshat) ./ bshatdenom).^2); %abs
        idx_calc = self.tfused_segs_sep_z;
    case {'HeteroEstVar', 'HeteroscedasticityEstVar'}

        if isempty(bshatvar)
            warning('bshatvar is empty, use bsvar')
            bshatvar = bsvar;
        end

        % bs = bs(idx_calc);
        % bshat = bshat(idx_calc);
        % bsstd = sqrt(bsvar(idx_calc) )+ self.dn_reg;%sqrt
        % bshatstd = sqrt(bshatvar(idx_calc) )+ self.dn_reg; %sqrt
        % if ~isempty(self.s_th)
        %     bshatstd(bshatstd < self.s_th) = 0;
        % end
        % loss = ([(real(bs) ./ bsstd - real(bshat) ./ bshatstd); ...
        %              (imag(bs) ./ bsstd - imag(bshat) ./ bshatstd)]); %abs
        % loss = (bs ./ bsstd) - (bshat ./ bshatstd);
        % loss = [real(loss); imag(loss)];

        bshatstd = sqrt(bshatvar) + self.dn_reg;
        % loss = sqrt([(real(bs) - real(bshat)) .^ 2 ./ bshatvar; (imag(bs) - imag(bshat)) .^ 2 ./ bshatvar]);
        % loss = sqrt(cat(3, (real(bs) - real(bshat)) .^ 2 ./ bshatvar, (imag(bs) - imag(bshat)) .^ 2 ./ bshatvar));
        % sep
        loss = cat(3, (real(bs) - real(bshat)) ./ bshatstd, (imag(bs) - imag(bshat)) ./ bshatstd);
        idx_calc = self.tfused_segs_sep_z;
        % complex
        % loss = abs(bs - bshat) .^ 2 ./ bshatvar;
        % idx_calc = self.tfused_segs;

    case {'Hetero', 'Heteroscedasticity'}
        % loss=(abs(bs(idx_calc)-bshat(idx_calc)).^2)./(bshatvar(idx_calc));
        % bs = bs(idx_calc);
        % bshat = bshat(idx_calc);
        % bshatvar=(bshatvar(idx_calc));
        % loss=sqrt([(real(bs)-real(bshat)).^2./bshatvar;(imag(bs)-imag(bshat)).^2./bshatvar]);

        % bsvar = (bsvar) + self.dn_reg;
        % loss = sqrt([(real(bs) - real(bshat)) .^ 2 ./ bsvar; (imag(bs) - imag(bshat)) .^ 2 ./ bsvar]);

        bsstd = sqrt(bsvar) + self.dn_reg;
        % loss = ([(real(bs) - real(bshat)) ./ bsstd; (imag(bs) - imag(bshat)) ./ bsstd]);
        loss = cat(3, (real(bs) - real(bshat)) ./ bsstd, (imag(bs) - imag(bshat)) ./ bsstd);
        idx_calc = self.tfused_segs_sep_z;
    case {'Leonenko'}
        % B2T=sqrt(self.Fs/2);
        % loss=B2T^2*(abs(bshat(idx_calc)-bs(idx_calc)).^2)./((W3.^2).*sxsysxyhat(idx_calc)./W2xW2yW2xy);%+1e-10
        % loss=B2T^2*(abs(bshat(idx_calc)-bs(idx_calc)).^2)./(bshatvar(idx_calc));
        % loss=(abs(bs(idx_calc)-bshat(idx_calc)).^2)./(bshatvar(idx_calc));
        % B2T=self.hos.NW/self.hos.N;
        % B2T=1;
        % loss=((B2T.^2)/2).*(abs(bs(idx_calc)-bshat(idx_calc)).^2)./(bshatvar(idx_calc));
        loss = (abs(bs(idx_calc) - bshat(idx_calc)) .^ 2) ./ (bshatvar(idx_calc));
        % loss=loss/length(idx_calc);
        % loss=(loss/(self.hos.N^2))/12;
        % loss=(loss/(self.hos.N));
    case {'mse'}
        loss = (abs(bshat(idx_calc) - bs(idx_calc)) .^ 2);
    case {'BR', 'AIC', 'BRNg'}

        switch lower(self.normalization)
            case {'haubrich'}
                loss = (abs(bs(idx_calc) - bshat(idx_calc)) .^ 2) ./ (bshatvar(idx_calc));
            case {'skewness'}
                loss = (abs(bs(idx_calc) - bshat(idx_calc)) .^ 2) ./ (bshatvar(idx_calc)); %+1e-10 (2*pi/self.Nf)*
            otherwise
                error('no such method')
        end

    case {'Rice'}
        idx_calc = idx_calc & bshat ~= 0 & bshat > (10 * eps);
        loss = abs((bs(idx_calc) - bshat(idx_calc)) ./ (bshat(idx_calc))) .^ 2; %any(isinf(loss))||any(isnan(loss))
        loss(isinf(loss)) = 0;
        % loss=abs((bs(idx_calc)-bshat(idx_calc))./(sqrt(bshatvar(idx_calc)))).^2;
        % loss=(abs(bs(idx_calc)-bshat(idx_calc)).^2)./(bshatvar(idx_calc));
    case {'MCE'}
        % idx_calc=idx_calc&~isnan(bshat)&bshat~=0&bshat>(10*eps);
        loss = log(abs(bshat(idx_calc)) .^ 2) + abs(bs(idx_calc)) .^ 2 ./ abs(bshat(idx_calc)) .^ 2;
    otherwise
        error('no such loss_type')
end

loss = loss(idx_calc);
loss(isnan(loss) | (isinf(loss) & loss > 0)) = inf;

if strcmpi(self.para_fit.optfun, 'fmincon')
    loss = sum(loss, "all");
end

% loss=(1./self.hos.N^2)*loss;
% loss=(1./self.hos.N)*loss;
% loss=(1./self.Nf)*loss;
loss = convert_to_double(loss);

end
