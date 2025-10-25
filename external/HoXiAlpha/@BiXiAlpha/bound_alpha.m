function [lb, ub] = bound_alpha(self)

if strcmpi(self.para_alpha.type, 'gaussian')
    %{
    s_min = 0;
    s_max = 1.1 * max(self.s, [], 'all');
    %     lb = [self.para_alpha.lb_mu, 2*self.df,...%mu, sigma
    %             0.01*self.para_xi.alpha,zeros(1,self.kmax-1)];%alpha
    if isempty(self.mufix)
        lb = [self.para_alpha.lb_mu, min(2 * self.df, 0.5), ... %mu, sigma
                  s_min .* ones(1, self.kmax)]; %alpha

        ub = [self.para_alpha.ub_mu, 2, ... %mu, sigma
                  1 * s_max .* ones(1, self.kmax)]; %alpha
    else
        lb = [self.mufix, min(2 * self.df, 0.5), ... %mu, sigma
                  s_min .* ones(1, self.kmax)]; %alpha

        ub = [self.mufix, 2, ... %mu, sigma
                  1 * s_max .* ones(1, self.kmax)]; %alpha
    end
    %}
    s_min = 0; %min(self.s);
    s_max = 2 * max(self.s, [], 'all');

    self.para_alpha.kernel.lb.h = s_min .* ones(self.kmax, 1);
    self.para_alpha.kernel.lb.mu = self.para_alpha.mu_init - 0.5;
    self.para_alpha.kernel.lb.sigma = 0.05 .* ones(self.kmax, 1); %0.75; % 0.001  0.5 0.5^2 min(2*self.df,0.1).^2;

    self.para_alpha.kernel.ub.h = s_max .* ones(self.kmax, 1);
    self.para_alpha.kernel.ub.mu = self.para_alpha.mu_init + 0.5;
    self.para_alpha.kernel.ub.sigma = 6 .* ones(self.kmax, 1);

    if ~isempty(self.para_alpha.kernel.fix)
        self.para_alpha.kernel.vlb(~isnan(self.para_alpha.kernel.vfix)) = self.para_alpha.kernel.vfix(~isnan(self.para_alpha.kernel.vfix));
        self.para_alpha.kernel.vub(~isnan(self.para_alpha.kernel.vfix)) = self.para_alpha.kernel.vfix(~isnan(self.para_alpha.kernel.vfix));
    end

    lb = self.para_alpha.kernel.vlb;
    ub = self.para_alpha.kernel.vub;
elseif strcmpi(self.para_alpha.type, 'tstudent') && strcmpi(self.peak_relation, 'harmonic')
    s_min = 0; %min(self.s);
    s_max = 1.1 * max(self.s, [], 'all');

    if self.para_alpha.no_zero_peak
        self.para_alpha.kernel.lb.h = s_min .* ones(self.kmax, 1);
    else
        self.para_alpha.kernel.lb.h = [0; s_min .* ones(self.kmax - 1, 1)];
    end

    % self.para_alpha.kernel.lb.h = s_min .* ones(self.kmax, 1);

    self.para_alpha.kernel.lb.mu = 4;
    self.para_alpha.kernel.lb.sigma = 0.3; %0.75; % 0.001  0.5 0.5^2 min(2*self.df,0.1).^2;
    self.para_alpha.kernel.lb.nu = 1; %1e3 10 1e-3 0.1
    self.para_alpha.kernel.lb.d = 2; %1e-3 1e-3 1e-3 1e-3
    self.para_alpha.kernel.lb.b = 1; %0

    if self.para_alpha.no_zero_peak
        self.para_alpha.kernel.ub.h = s_max .* ones(self.kmax, 1);
    else
        self.para_alpha.kernel.ub.h = [0; s_max .* ones(self.kmax - 1, 1)];
    end

    % self.para_alpha.kernel.ub.h = s_max .* ones(self.kmax, 1);

    self.para_alpha.kernel.ub.mu = 20;
    % self.para_alpha.kernel.ub.sigma = 10; %2.5 2;
    self.para_alpha.kernel.ub.sigma = 2.75;
    self.para_alpha.kernel.ub.nu = 20;
    self.para_alpha.kernel.ub.d = 4;
    self.para_alpha.kernel.ub.b = 1; %1e5

    self.para_alpha.kernel.lb.h(1) = 0;
    self.para_alpha.kernel.ub.h(1) = 0;

    if ~isempty(self.mufix)
        self.para_alpha.kernel.lb.mu = self.mufix;
        self.para_alpha.kernel.ub.mu = self.mufix;
    elseif isempty(self.mufix) && strcmpi(self.order, '1')
        % self.para_alpha.kernel.lb.sigma=self.para_alpha.kernel.para.sigma;
        % self.para_alpha.kernel.ub.sigma=self.para_alpha.kernel.para.sigma;
        % self.para_alpha.kernel.lb.nu=self.para_alpha.kernel.para.nu;
        % self.para_alpha.kernel.ub.nu=self.para_alpha.kernel.para.nu;
    end

    % if strcmpi(self.order,'1')
    %     self.para_alpha.kernel.lb.mu=self.para_alpha.kernel.para.mu;
    %     self.para_alpha.kernel.ub.mu=self.para_alpha.kernel.para.mu;
    % end
    if ~isempty(self.para_alpha.kernel.fix)
        self.para_alpha.kernel.vlb(~isnan(self.para_alpha.kernel.vfix)) = self.para_alpha.kernel.vfix(~isnan(self.para_alpha.kernel.vfix));
        self.para_alpha.kernel.vub(~isnan(self.para_alpha.kernel.vfix)) = self.para_alpha.kernel.vfix(~isnan(self.para_alpha.kernel.vfix));
    end

    lb = self.para_alpha.kernel.vlb;
    ub = self.para_alpha.kernel.vub;
elseif strcmpi(self.para_alpha.type, 'tstudent') && strcmpi(self.peak_relation, 'free')
    s_min = 0; %min(self.s);
    s_max = 1.1 * max(self.s, [], 'all');

    if self.para_alpha.no_zero_peak
        self.para_alpha.kernel.lb.h = s_min .* ones(self.kmax, 1);
    else
        self.para_alpha.kernel.lb.h = [0; s_min .* ones(self.kmax - 1, 1)];
    end

    % self.para_alpha.kernel.lb.h = s_min .* ones(self.kmax, 1);

    % self.para_alpha.kernel.lb.mu = 4.* ones(self.kmax, 1);
    self.para_alpha.kernel.lb.mu = self.para_alpha.mu_init - 0.5;
    self.para_alpha.kernel.lb.sigma = 0.3 .* ones(self.kmax, 1); %0.75; % 0.001  0.5 0.5^2 min(2*self.df,0.1).^2;
    self.para_alpha.kernel.lb.nu = 1 .* ones(self.kmax, 1); %1e3 10 1e-3 0.1
    self.para_alpha.kernel.lb.d = 2 .* ones(self.kmax, 1); %1e-3 1e-3 1e-3 1e-3
    % self.para_alpha.kernel.lb.b = 0.5 .* ones(self.kmax, 1); %0
    self.para_alpha.kernel.lb.b = 1 .* ones(self.kmax, 1); %0

    if self.para_alpha.no_zero_peak
        self.para_alpha.kernel.ub.h = s_max .* ones(self.kmax, 1);
    else
        self.para_alpha.kernel.ub.h = [0; s_max .* ones(self.kmax - 1, 1)];
    end

    % self.para_alpha.kernel.ub.h = s_max .* ones(self.kmax, 1);

    % self.para_alpha.kernel.ub.mu = 20 .* ones(self.kmax, 1);
    self.para_alpha.kernel.ub.mu = self.para_alpha.mu_init + 0.5; % self.para_alpha.kernel.ub.sigma = 10; %2.5 2;
    self.para_alpha.kernel.ub.sigma = 2.75 .* ones(self.kmax, 1);
    self.para_alpha.kernel.ub.nu = 20 .* ones(self.kmax, 1);
    self.para_alpha.kernel.ub.d = 4 .* ones(self.kmax, 1);
    % self.para_alpha.kernel.ub.b = 2 .* ones(self.kmax, 1); %1e5
    self.para_alpha.kernel.ub.b = 1 .* ones(self.kmax, 1); %1e5

    if ~isempty(self.mufix)
        self.para_alpha.kernel.lb.mu = self.mufix;
        self.para_alpha.kernel.ub.mu = self.mufix;
    elseif isempty(self.mufix) && strcmpi(self.order, '1')
        % self.para_alpha.kernel.lb.sigma=self.para_alpha.kernel.para.sigma;
        % self.para_alpha.kernel.ub.sigma=self.para_alpha.kernel.para.sigma;
        % self.para_alpha.kernel.lb.nu=self.para_alpha.kernel.para.nu;
        % self.para_alpha.kernel.ub.nu=self.para_alpha.kernel.para.nu;
    end

    if ~isempty(self.para_alpha.kernel.fix)
        self.para_alpha.kernel.vlb(~isnan(self.para_alpha.kernel.vfix)) = self.para_alpha.kernel.vfix(~isnan(self.para_alpha.kernel.vfix));
        self.para_alpha.kernel.vub(~isnan(self.para_alpha.kernel.vfix)) = self.para_alpha.kernel.vfix(~isnan(self.para_alpha.kernel.vfix));
    end

    lb = self.para_alpha.kernel.vlb;
    ub = self.para_alpha.kernel.vub;
end

end

% find_max_around_xt(self.f,self.s,self.para_alpha,mu)
