function show(self)

% if ~isempty(self.bs)
%     % tiledlayout(3, 4)
%     % subplot(3,4)
% end

% set(gcf, 'Position', get(0, 'Screensize'));
pause(0.00001);
% frame_h = get(handle(gcf),'JavaFrame');
% set(frame_h,'Maximized',1);
hFig = gcf;
clf
hFig.WindowState = 'maximized';

if ~isempty(self.bs)
    % nexttile(1, [1, 4])
    subplot(3,1,1)
end

% plot real data, deal with multiple segments
cm = colormap('lines');
hold on;
% plot(self.f, 10 * log10(self.s), '-', 'Color', cm(1, :), 'LineWidth', 0.5, 'HandleVisibility', 'off');
plot(NaN, 'DisplayName', 'real data', 'Color', cm(1, :));
plot(self.f, 10 * mean(log10(self.s), 2), '-x', 'Color', cm(2, :), 'LineWidth', 2);
plot(self.f, 10 * log10(self.shat), '-o', 'Color', cm(3, :), 'LineWidth', 2);

[~, scomp] = pred_s(self.f, self.ks, self.kh, self.para_xi, self.para_alpha);

if strcmpi(self.para_xi.type, 'tstudent+peak')
    scomp(:, 1) = scomp(:, 1) + scomp(:, 2);
    scomp(:, 2) = [];
end

smin = min(10 * [mean(log10(self.s), 2); log10(self.shat)], [], 'all');
smax = max(10 * [mean(log10(self.s), 2); log10(self.shat)], [], 'all');
% scomp(10 * log10(scomp) < smin - 2) = nan;

% plot(self.f, 10 * log10(scomp), '-', 'LineWidth', 2); hold on
for icomp = 1:size(scomp, 2)
    % if all(isnan(scomp(:,icomp)))
    %     continue;
    % end
    plot(self.f, 10 * log10(scomp(:, icomp)), '-', 'Color', cm(3 + icomp, :), 'LineWidth', 2); hold on
end
if strcmpi(self.peak_relation, 'harmonic')
    mu = get_harmonic_mu_sigma(self.para_alpha.kernel.para.mu, 1, self.ks, self.kh, self.para_alpha.no_zero_peak);
elseif strcmpi(self.peak_relation, 'free')
    mu = self.para_alpha.kernel.para.mu(:)';
end
xline(mu, ':k', 'Alpha', 0.6, 'LineWidth', 5, 'HandleVisibility', 'off')

xlim([0, max(self.f)])
% ylim([smin - 5, smax + 5])
% ylim([max(smin - 5, min(smax + 5, - 10)), smax + 5])
% ylims=ylim();
ylim([smin - 2, smax+2])
xticks([0:1:max(self.f)])

if self.ks > 0
    lagalpha = [append('1/', string(num2str((self.ks + 1:-1:2)')), ' harmonic'); "alpha peak"; append(string(num2str((2:self.kh)')), '/1 harmonic')];
else
    lagalpha = ["alpha peak"; append(string(num2str((2:self.kh)')), '/1 harmonic')];
end

if ~self.para_alpha.no_zero_peak
    lagalpha = ['0 peak'; lagalpha]; %\delta f peak
end

lagempty = repmat('', 1, self.ks + self.kh);
ws = warning('off');
legend(["real data", "mean of data", "fitted data", "\xi process", lagalpha(1:self.ks + self.kh)', lagempty])
warning(ws);
title('spectrum')

if isempty(self.bs)
    return;
end

%% third order part
bc = mean(self.bc, 3);
bchat = mean(self.bchat, 3);
bc_xi = mean(self.bchat_component('xi'), 3);
bc_alpha = mean(self.bchat_component('alpha'), 3);

dtick = 5;
range = [0:dtick:max(self.f)];

% nexttile(5, [1, 1])
subplot(3,4,5)
surfbc(self.fx, self.fy, abs(bc));
xticks(range)
yticks(range)
hold on; plot_bsgrid(range)
title(['orignal bicoherence [', self.normalization, ']'])

% nexttile(9, [1, 1])
subplot(3,4,9)
surfbc(self.fx, self.fy, abs(bchat));
xticks(range)
yticks(range)
hold on; plot_bsgrid(range)
title(['fitted bicoherence [', self.normalization, ']'])

% warning('plotting bispectrum')
% nexttile(6, [1, 1])
subplot(3,4,6)
surfbs(self.fx, self.fy, real(bc));
xticks(range)
yticks(range)
hold on; plot_bsgrid(range)
title('real part of orignal bicoherence')

% nexttile(10, [1, 1])
subplot(3,4,10)
surfbs(self.fx, self.fy, real(bchat));
% clim([min(real(bc(~isinf(bc)))),max(real(bc(~isinf(bc))))])
xticks(range)
yticks(range)
hold on; plot_bsgrid(range)
title('real part of fitted bicoherence')

% nexttile(7, [1, 1])
subplot(3,4,7)
surfbs(self.fx, self.fy, imag(bc));
xticks(range)
yticks(range)
hold on; plot_bsgrid(range)
title('imaginary part of orignal bicoherence')

% nexttile(11, [1, 1])
subplot(3,4,11)
surfbs(self.fx, self.fy, imag(bchat));
% clim([min(imag(bc(~isinf(bc)))),max(imag(bc(~isinf(bc))))])
xticks(range)
yticks(range)
hold on; plot_bsgrid(range)
title('imaginary part of fitted bicoherence')

%%
%%{
% nexttile(8, [1, 1])
subplot(3,4,8)
surfbs(self.fx, self.fy, abs(bc_xi));
xticks(range)
yticks(range)
hold on; plot_bsgrid(range)
title('bicoherence of \xi process')

% nexttile(12, [1, 1])
subplot(3,4,12)
surfbc(self.fx, self.fy, abs(bc_alpha));
% clim([min(imag(bc(~isinf(bc)))),max(imag(bc(~isinf(bc))))])
xticks(range)
yticks(range)
hold on; plot_bsgrid(range)
title('bicoherence of \alpha process')
%%}
%{
nexttile(8, [1, 1])
surfbs(self.fx, self.fy, abs(self.bchat_component('xi+xi')));
xticks(range)
yticks(range)
hold on; plot_bsgrid(range)
title('bicoherence of \xi process')

nexttile(12, [1, 1])
surfbs(self.fx, self.fy, abs(self.bchat_component('alpha+alpha')));
% clim([min(imag(bc(~isinf(bc)))),max(imag(bc(~isinf(bc))))])
xticks(range)
yticks(range)
hold on; plot_bsgrid(range)
title('bicoherence of \alpha process')
%}
%{
nexttile(8, [1, 1])
surfbs(self.fx, self.fy, log10(abs(self.bshat_component('xi'))));
xticks(range)
yticks(range)
hold on; plot_bsgrid(range)
title('bicoherence of \xi process')
% z=zlim;

nexttile(12, [1, 1])
surfbs(self.fx, self.fy, log10(abs(self.bshat_component('alpha'))));
% clim([min(imag(bc(~isinf(bc)))),max(imag(bc(~isinf(bc))))])
xticks(range)
yticks(range)

hold on; plot_bsgrid(range)
% zlim(z)
title('bicoherence of \alpha process')
%}
%{
nexttile(8, [1, 1])
% surfbs(self.fx, self.fy, abs(self.bshat_component('xi+xi')));
surfbs(self.fx, self.fy, log10(abs(self.bshat)));

xticks(range)
yticks(range)
hold on; plot_bsgrid(range)
title('bispectrum')

nexttile(12, [1, 1])
surfbs(self.fx, self.fy, log10(abs(self.sxsysxy)));
% clim([min(imag(bc(~isinf(bc)))),max(imag(bc(~isinf(bc))))])
xticks(range)
yticks(range)
hold on; plot_bsgrid(range)
title('sxsysxy')
%}

%{
nexttile(8, [1, 1])
% surfbs(self.fx, self.fy, abs(self.bshat_component('xi+xi')));
surfbs(self.fx, self.fy, log10(abs(self.bs)));

xticks(range)
yticks(range)
hold on; plot_bsgrid(range)
title('bispectrum')

nexttile(12, [1, 1])
surfbs(self.fx, self.fy, log10(abs(self.bshat)));
% clim([min(imag(bc(~isinf(bc)))),max(imag(bc(~isinf(bc))))])
xticks(range)
yticks(range)
hold on; plot_bsgrid(range)
title('bispectrum fitted')
%}

if ~isempty(self.name)
    sgtitle(replace(append(self.name, ['-', self.hos_components]), '_', '-'))
end

%{
figure
tiledlayout(3, 2)
bs = mean(self.bs, 3);
bshat = mean(self.bshat, 3);
bs_xi = mean(self.bshat_component('xi'), 3);
bs_alpha = mean(self.bshat_component('alpha'), 3);

nexttile(1, [1, 1])
% surfbs(self.fx, self.fy, abs(self.bshat_component('xi+xi')));
surfbs(self.fx, self.fy, log10(abs(bs)));
xticks(range)
yticks(range)
hold on; plot_bsgrid(range)
title('bispectrum')
z = zlim;
c = clim;

nexttile(2, [1, 1])
surfbs(self.fx, self.fy, log10(abs(bshat)));
% clim([min(imag(bc(~isinf(bc)))),max(imag(bc(~isinf(bc))))])
xticks(range)
yticks(range)
hold on; plot_bsgrid(range)
title('bispectrum fitted')
zlim(z)
clim(c)

nexttile(3, [1, 1])
surfbs(self.fx, self.fy, log10(abs(bs_xi)));
% clim([min(imag(bc(~isinf(bc)))),max(imag(bc(~isinf(bc))))])
xticks(range)
yticks(range)
hold on; plot_bsgrid(range)
title('bispectrum fitted- \xi process')
zlim(z)
clim(c)

nexttile(4, [1, 1])
surfbs(self.fx, self.fy, log10(abs(bs_alpha)));
% clim([min(imag(bc(~isinf(bc)))),max(imag(bc(~isinf(bc))))])
xticks(range)
yticks(range)
hold on; plot_bsgrid(range)
title('bispectrum fitted- \alpha process')
zlim(z)
clim(c)

nexttile(5, [1, 1])
surfbs(self.fx, self.fy, log10(abs(mean(self.bsdenom, 3))));
% clim([min(imag(bc(~isinf(bc)))),max(imag(bc(~isinf(bc))))])
xticks(range)
yticks(range)
hold on; plot_bsgrid(range)
title('denominator for bicoherence')
zlim(z)
clim(c)

nexttile(6, [1, 1])
surfbs(self.fx, self.fy, log10(abs(self.bshatdenom)));
% clim([min(imag(bc(~isinf(bc)))),max(imag(bc(~isinf(bc))))])
xticks(range)
yticks(range)
hold on; plot_bsgrid(range)
title('denominator fitted for bicoherence')
zlim(z)
clim(c)

%}

end
