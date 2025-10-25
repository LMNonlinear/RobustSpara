function [smax,fmax]=findmax_spectrum(s,f,range)
idx=find(f>=range(1)&f<=range(2));
f=f(idx,:);
s=s(idx,:);
% s_mirror=s;
% s=[s_mirror;s;s_mirror];
% [y(4:-1:2) y y(end-1:-1:end-3)]
% findpeaks(s)
[smax,idx_max]=max(s);
fmax=f(idx_max);


end


