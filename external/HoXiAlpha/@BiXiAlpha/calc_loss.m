function [loss, p, q] = calc_loss(self, loss_type, s, shat, svar, bs, bshat, bsvar, bshatvar)

if isempty(loss_type)
    loss_type = self.loos_type;
end

if strcmpi(self.para_fit.optfun, 'fmincon')
    loss = self.p * calc_loss_s(self, loss_type, s, shat, svar) + self.q * calc_loss_bs(self, loss_type, bs, bshat, bsvar, bshatvar);

else
    %{
    % loss = [self.p .* calc_loss_s(self, loss_type, s, shat, svar); self.q .* calc_loss_bs(self, loss_type, bs, bshat, bsvar, bshatvar)];
    loss_s = calc_loss_s(self, loss_type, s, shat, svar);
    loss_bs = calc_loss_bs(self, loss_type, bs, bshat, bsvar, bshatvar);
    Ns = numel(loss_s);
    Nbs = numel(loss_bs);
    p = self.p .* (Nbs / Ns);
    q = self.q;
    loss = [p .* loss_s; q .* loss_bs];
    %}
    % loss = [self.p .* calc_loss_s(self, loss_type, s, shat, svar); self.q .* calc_loss_bs(self, loss_type, bs, bshat, bsvar, bshatvar)];
    loss_s = calc_loss_s(self, loss_type, s, shat, svar);
    loss_bs = calc_loss_bs(self, loss_type, bs, bshat, bsvar, bshatvar);

    % loss = [self.p .* loss_s; self.q .* loss_bs];
    loss = [loss_s; loss_bs];

end

loss = convert_to_double(loss);
end
