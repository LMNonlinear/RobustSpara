function x=vec_para(self,part,ispad2max)
if nargin<3 ||isempty(ispad2max)
    ispad2max=true;
end

switch part
    case {'k'}
        x=[self.ks;self.kh];
    case {'xi'}
        % switch self.para_xi.type
        %     case{'tstudent'}
        %         x=self.para_xi.kernel.vpara();
        %     otherwise
        %         error('no such method')
        % end
        x=self.para_xi.kernel.vpara();
    case{'alpha'}
        x=self.para_alpha.kernel.vpara();
    case{'bs'}
        x=self.para_bs.kernel.vpara();
    otherwise
        error('no such part')
end
end