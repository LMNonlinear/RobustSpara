function surfbc(fx, fy, bc)

p = surf(fx, fy, bc);
view(0, 90)
xlim([min(fx, [], 'all'), max(fx, [], 'all')])
ylim([min(fy, [], 'all'), max(fy, [], 'all')])
set(p, 'edgecolor', 'none')
% shading interp
grid off
ax = gca;
% axis equal
% colormap(turbo)% colormap(flip(hot))%hot
colormap(ax, flip(hot))
colorbar
title('bicoherence')
box off
clim([0.05, 0.6])

end
