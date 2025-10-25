function init_xi(self)
%% fit the trend
% self.trend=fit_trend(self);
%% init xi para
init_para_xi(self);
bound_xi(self);
% %% get the better init para from the trend
% if strcmpi(self.para_xi.type,'tstudent') || strcmpi(self.para_xi.type,'tstudent+peak')
%     range=0+[-self.width_candidate,self.width_candidate];
%     self.para_xi.kernel.para.h     = findmax_spectrum(self.s,self.f,range);
%     self.para_xi.alpha = self.para_xi.kernel.para.h./10;
% 
%     % find frequency in half maximum of the trend to get the parameter B of the xi model
%     % idx=find_closest(self.trend,max(self.trend)/2);% [idx]=find_closest(10.^trend,self.para_xi.kernel.para.h/2);
%     % self.para_xi.fhm=self.f(idx);% frequency in half maximum
%     % self.para_xi.kernel.para.sigma=self.para_xi.fhm/(sqrt(2^(1/3.2)-1));
%     self.para_xi.kernel.para.sigma=50;
% 
% elseif strcmpi(self.para_xi.type,'linear')
%     self.para_xi.w = polyfit(self.f, self.s, 1); 
% else
%     error('no such method')
% end

fit_xi(self);
end





% plot(self.f,self.s);hold on;
% plot(self.f,self.trend);
% plot(self.f,studentt_waveform(self.f,self.para_xi.kernel.para.h,self.para_xi.kernel.para.sigma,self.para_xi.kernel.para.nu,self.para_xi.kernel.para.mu));hold on









