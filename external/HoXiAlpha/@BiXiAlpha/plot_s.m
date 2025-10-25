function plot_s_shat(self)

pause(0.00001);
% frame_h = get(handle(gcf),'JavaFrame');
% set(frame_h,'Maximized',1);
hFig = gcf;
% hFig.WindowState = 'maximized';

FontSize=6;
LineWidth=1;
MarkerSize=2;
set(0,'defaultAxesFontName', 'arial')
set(0,'defaultTextFontName', 'arial')


p = get(0, "MonitorPositions");
% hFig.Position = p(2, :);
% hFig.Position= [200 200 600 300];
Position=p(2,:);
Position=[Position(1)+200,Position(2)+200,270,150];
hFig.Position=Position;






%%

% plot real data, deal with multiple segments
% cm = colormap('lines');
% cm=flip(plt.Accent);
cm=(plt.Set1);
hold on;
% plot(self.f, 10 * log10(self.s), '-', 'Color', cm(1, :), 'LineWidth', 0.5, 'HandleVisibility', 'off');
plot(NaN, 'DisplayName', 'real data', 'Color', cm(1, :));
% plot(self.f, 10 * log10(mean(self.s, 2)), '-o', 'Color', cm(2, :), 'LineWidth', LineWidth,'MarkerSize',MarkerSize);
% plot(self.f, 10 * log10(self.shat), '-x', 'Color', cm(3, :), 'LineWidth', LineWidth,'MarkerSize',MarkerSize);
plot(self.f, 10 * log10(mean(self.s, 2)), 'Color', cm(2, :), 'LineWidth', LineWidth,'MarkerSize',MarkerSize);

smin = min(10 * [log10(mean(self.s, 2)); log10(self.shat)], [], 'all');
smax = max(10 * [mean(log10(self.s), 2); log10(self.shat)], [], 'all');


xlim([0, max(self.f)])
ylim([smin - 2, smax + 2])
% ylim([max(smin - 5, min(smax + 5, - 10)), smax + 5])
xticks([0:5:max(self.f)])
xlabel('Frequency (Hz)')


ylabel('Power Spectrum density (db)')


hA=gca;
hA.XAxis.FontSize=FontSize;
hA.YAxis.FontSize=FontSize;
hA.YRuler.TickLabelGapOffset=5;
xtickangle(45);

lgd=legend(["", "empirical spectra"],'FontSize',FontSize);
lgd.FontName='Arial';

drawnow; pause(0.05);  % this innocent line prevents the Matlab hang

end