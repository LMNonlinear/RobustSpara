function [lb,ub]=bound_k(self)

% lb=self.kh;
% ub=self.kh;
% if isfield(self,'para_alpha')&& isfield(self.para_alpha,'alpha') && ~isempty(self.para_alpha.kernel.para.mu)
%     warning('test kh==5')
%     lb=5;
%     ub=ceil(max(self.f)./4);% ub=ceil(max(self.f)./self.para_alpha.kernel.para.mu);
% else
%     lb=0;
%     ub=8;
% end
% warning('test kh==5')
% lb=0;
% ub=ceil(max(self.f)./4);
% lb=ceil(max(self.f)./4);
% ub=ceil(max(self.f)./4);
if strcmpi(self.ktype, 'fix')
    lb=[self.ks;self.kh];
    ub=[self.ks;self.kh];
elseif strcmpi(self.ktype, 'free')
    lb=[0;min(ceil(max(self.f)/self.para_alpha.kernel.para.mu),self.khmax)];
    % ub=[ceil(self.para_alpha.ub_mu./self.para_alpha.lb_mu),ceil(max(self.f)./self.para_alpha.lb_mu)];
    ub=[self.ksmax;self.khmax];
elseif strcmpi(self.ktype, 'max')
    lb=[self.ksmax;self.khmax];
    ub=[self.ksmax;self.khmax];
end
end







