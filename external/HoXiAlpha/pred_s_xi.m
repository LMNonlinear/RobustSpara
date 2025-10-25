function [shat, scomp] = pred_s_xi(f, para_xi)
%     if strcmpi(type_xi,'tstudent') || strcmpi(type_xi,'tstudent+peak')
%         % peak is very sharp so only consider the tail of tstudent
%         self.alpha(1:end) = self.alpha(1:end) - studentt_waveform(mu1(1:end), self.para_xi.kernel.para.h, self.para_xi.kernel.para.sigma, self.para_xi.kernel.para.nu, self.para_xi.kernel.para.mu);
%     end

% if nargin < 3 || isempty(para_bs)

%     if strcmpi(para_xi.type, 'tstudent')
%         scomp = studentt_waveform(f, para_xi.kernel.para.h, para_xi.kernel.para.sigma, para_xi.kernel.para.nu, para_xi.kernel.para.mu);
%     elseif strcmpi(para_xi.type, 'tstudent+peak')
%         scomp(:, 1) = studentt_waveform(f, para_xi.kernel.para.h, para_xi.kernel.para.sigma, para_xi.kernel.para.nu, para_xi.kernel.para.mu);
%         scomp(:, 2) = gaussian_kernel(f, para_xi.kernel.para.mu, para_xi.kernel.para.sigma) * para_xi.alpha;
%     % elseif strcmpi(para_xi.type, 'exponential')
%     %     scomp(:, 1) = exponential_kernel(f, para_xi.kernel.para.h, self.sigma, para_xi.kernel.para.nu);
%     elseif strcmpi(para_xi.type,'gaussian')
%         scomp = para_xi.kernel.eval(f);
%     elseif strcmpi(self.para_xi.type,'linear')
%         scomp = polyval(self.para_xi.w, f);
%     else
%         error('no such method')
%     end

% else % for the pred_bs with different B(bandwidth) of xi

%     if strcmpi(para_xi.type, 'tstudent')
%         scomp = studentt_waveform(f, para_xi.kernel.para.h, para_bs.xi_B, para_xi.kernel.para.nu, para_xi.kernel.para.mu);

%     elseif strcmpi(para_xi.type, 'tstudent+peak')
%         scomp(:, 1) = studentt_waveform(f, para_xi.kernel.para.h, para_bs.xi_B, para_xi.kernel.para.nu, para_xi.kernel.para.mu);
%         scomp(:, 2) = gaussian_kernel(f, para_xi.kernel.para.mu, para_bs.xi_sigma) * para_xi.alpha;
%     elseif strcmpi(para_xi.type,'gaussian')
%         para_xi.kernel.para.sigma=para_bs.xi_sigma;
%         scomp = para_xi.kernel.eval(f);
%     elseif strcmpi(self.para_xi.type,'linear')
%         scomp = polyval(self.para_xi.w, f);
%     else
%         error('no such method')
%     end

% end
if nargout > 1
    [shat, scomp] = para_xi.kernel.eval(f);
else
    [shat] = para_xi.kernel.eval(f);
    scomp = [];
end

% shat = sum(scomp, 2);

end
