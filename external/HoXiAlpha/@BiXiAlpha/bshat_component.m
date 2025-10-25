function bshat = bshat_component(self, component, para_bs)

if nargin < 3 || isempty(para_bs)
    para_bs = self.para_bs;
end

switch component
    case {'xi','xi_xi'}
        mask = false(size(self.para_bs.kernel.para.h));
        mask(1, :, :) = true;
        mask(:, 1, :) = true;
        para_bs.kernel.para.h(~mask) = 0;
        % para_bs.kernel.para.h(:,:,2)(~mask)=0;
        % if para_bs.additive
        %     para_bs.betaxy(2:end)=0;
        %     para_bs.gammaxy(2:end)=0;
        % end
    case {'alpha','alpha_alpha'}
        mask = true(size(self.para_bs.kernel.para.h));
        mask(1, 1, :) = false;
        para_bs.kernel.para.h(~mask) = 0;
        % para_bs.kernel.para.h(:,:,2)(~mask)=0;
        % if para_bs.additive
        %     para_bs.betaxy(1)=0;
        %     para_bs.gammaxy(1)=0;
        % end
    case {'xi+xi'}
        mask = false(size(self.para_bs.kernel.para.h));
        mask(1, 1, :) = true;
        para_bs.kernel.para.h(~mask) = 0;
        % para_bs.kernel.para.h(~mask)=0;
        % if para_bs.additive
        %     para_bs.betaxy(2:end)=0;
        %     para_bs.gammaxy(2:end)=0;
        % end
    case {'alpha+alpha'}
        mask = true(size(self.para_bs.kernel.para.h));
        mask(1, :, :) = false;
        mask(:, 1, :) = false;
        para_bs.kernel.para.h(~mask) = 0;
        % para_bs.kernel.para.h(:,:,2)(~mask)=0;
        % if para_bs.additive
        %     para_bs.betaxy(1)=0;
        %     para_bs.gammaxy(1)=0;
        % end
    case {'xi+alpha'}
        mask = false(size(self.para_bs.kernel.para.h));
        mask(1, :, :) = true;
        mask(:, 1, :) = true;
        mask(1, 1, :) = false;
        para_bs.kernel.para.h(~mask) = 0;
        % para_bs.kernel.para.h(:,:,2)(~mask)=0;
        % if para_bs.additive
        %     para_bs.betaxy(2:end)=0;
        %     para_bs.gammaxy(2:end)=0;
        % end
    otherwise
        error('the component is not exsit, only xi and alpha');
end

bshat = pred_bs(self.f, para_bs, self.fxfy);

end

% para_bs.betaxy(2:end)=0;
% para_bs.gammaxy(2:end)=0;
