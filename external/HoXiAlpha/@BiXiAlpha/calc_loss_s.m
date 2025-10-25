function loss = calc_loss_s(self, loss_type, s, shat, svar, shatvar)

if isempty(loss_type)
    loss_type = self.loos_type;
end

if nargin < 6 || isempty(shatvar)
    shatvar = self.shatvar;
end

if nargin < 5 || isempty(svar)
    svar = self.svar;
end

switch loss_type
    case {'bicoherence'}
        % loss = sqrt((log(shat) - log(s)) .^ 2 ./ svar);
        % loss = (abs(log(shat) - log(s)) .^ 2 ./ svar);
        loss = ((log(shat) - log(s)) ./ sqrt(svar));
    case {'Hetero', 'Heteroscedasticity', 'HeteroEstVar', 'HeteroscedasticityEstVar'}
        shat(shat < eps) = eps;
        % loss = sqrt((log(shat) - log(s)) .^ 2 ./ shatvar); % log transform make homosdestic

        % sqrt
        loss = ((log(shat) - log(s)) ./ sqrt(svar)); % log transform make homosdestic
        % sqr
        % loss = abs(log(shat) - log(s)).^2 ./ (svar); % log transform make homosdestic

        % case {'HeteroEst', 'HeteroEstVar', 'HeteroscedasticityEstVar'}
        %
        %     if isempty(shatvar)
        %         warning('shatvar is empty, use homoskedasticity')
        %         shatvar = ones(size(s));
        %     end
        %
        %     shatvar = shatvar + self.dn_reg;
        %     % loss = sqrt((log(shat) - log(s)) .^ 2 ./ shatvar); % log transform make homosdestic
        %     loss = ((log(shat) - log(s)) ./ sqrt(shatvar)); % log transform make homosdestic
    case 'Leonenko'
        % B1T=self.Fs/2;
        % loss=((shat-s)./shat).^2;
        % loss=B1T*loss;
        % idx_calc=~isnan(shat)&shat~=0;
        % loss=((shat(idx_calc)-s(idx_calc))./shat(idx_calc)).^2;
        % loss(isinf(loss))=0;
        % B1T=self.hos.NW/self.hos.N;
        % B1T=1;
        % idx_calc=~isnan(shat)&shat~=0;
        % loss=B1T.*((shat(idx_calc)-s(idx_calc))./shat(idx_calc)).^2;
        % loss(isinf(loss))=0;

        % idx_calc=~isnan(shat)&shat~=0;
        % loss=((shat(idx_calc)-s(idx_calc))./shat(idx_calc)).^2;
        % loss(isinf(loss))=0;
        % loss=loss/length(idx_calc);
        % loss=(loss/self.hos.N)/2;
        loss = ((shat - s) ./ shat) .^ 2;
    case 'mse'
        loss = (shat - s) .^ 2;
    case {'BR', 'AIC', 'MCE'} %,'BRabs'
        idx_calc = ~isnan(shat) & shat ~= 0;
        loss = log(shat(idx_calc)) + s(idx_calc) ./ shat(idx_calc); %((2*pi)^(-2)/self.Nf)*
        % loss=shat(idx_calc).*exp(s(idx_calc)./shat(idx_calc));
    case {'BRNg', 'BrillingerNongaussian'} %,'BRabs'
        idx_calc = ~isnan(shat) & shat ~= 0;
        % loss=log(sqrt(shat(idx_calc)))+0.5*sqrt((s(idx_calc)-shat(idx_calc)).^2)./shat(idx_calc);%((2*pi)^(-2)/self.Nf)*
        % loss=log(sqrt(shat(idx_calc)))+0.5*((s(idx_calc)-shat(idx_calc)).^2)./shat(idx_calc);
        % loss=0.5*log((shat(idx_calc)))+0.5*s(idx_calc)./shat(idx_calc);
        % loss=0.5*log(shat(idx_calc))+0.5*((s(idx_calc)-shat(idx_calc)).^2)./shat(idx_calc);
        loss = ((s(idx_calc) - shat(idx_calc)) .^ 2) ./ shat(idx_calc);
    case {'Rice'}
        idx_calc = ~isnan(shat) & shat ~= 0;
        loss = ((s(idx_calc) - shat(idx_calc)) ./ shat(idx_calc)) .^ 2;
    otherwise
        error('no such loss_type')
end

loss = loss(:);
loss(isnan(loss) | isinf(loss)) = inf;

if strcmpi(self.para_fit.optfun, 'fmincon')
    loss = sum(loss, "all");
end

% loss=sum(loss,"all");
% loss=(1./self.hos.N)*loss;
% loss=(1./self.Nf)*loss;
loss = convert_to_double(loss);

end
