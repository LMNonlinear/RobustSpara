function bc = calc_bc(self, bs, denominator)

switch self.normalization
    case 'skewness'
        % skewness function [Collis et al. (1998)]
        % bsvar=sxsysxy
        % bchat = sqrt(abs(self.bshat) .^ 2 ./ (self.shat .* self.shat.' .* self.sxyhat)); %
        bc = sqrt(abs(bs) .^ 2 ./ (denominator .^ 2 + self.dn_reg)); %+self.dn_reg
    case {'haubrich'}
        % Haubrich normalization which has no upper bound
        % [see Birkelund et al. (2001) and Haubrich (1965)]
        % bchat = self.bshat ./ (sqrt(self.shat .* self.shat.' .* self.sxyhat));
        % bsvar=sxsysxy
        bc = bs ./ ((denominator) + self.dn_reg); %+self.dn_reg

        %case {'shahbazi'}
        %    % Shahbazi et al. (2014)
        %   bchat = self.bshat ./ (self.shat .* self.shat.' .* self.sxyhat);
    case {'hagihira', 'hag'}
        bc = bs ./ ((denominator) + self.dn_reg);
    otherwise
        error('no such normalization')
end

end
