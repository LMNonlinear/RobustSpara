function plot_bc(self,comp,op,cmpath)
% if nargin<2||isempty(type)
%     type='bc';
% end

if nargin<4||isempty(cmpath)
    cmpath = [];
end
hFig = gcf;
% hFig.WindowState = 'maximized';

FontSize=6;
LineWidth=1;
p = get(0, "MonitorPositions");
% hFig.Position = p(2, :);
% hFig.Position= [200 200 600 300];
Position=p(2,:);
Position=[Position(1)+200,Position(2)+200,160,160];
% Position=[Position(1)+200,Position(2)+200,400,400];

hFig.Position=Position;


% plot_bc(self)
% plot_bchat(self)
% surfbc(self.fx,self.fy,abs(self.bc));
bc=op(self.(comp));
bc(self.fx+self.fy>max(self.f))=NaN;
surfbc(self.fx,self.fy,bc);


if sum(bc>0,"all")&&sum(bc<0,"all")
    clim([-0.6, 0.6]);
    % clm.Ticks = -0.6:0.3:0.6;
    % AxesH = axes('CLim', [-0.6, 0.6]);
    % clm = colorbar('peer', AxesH, 'h', ...
    %     'XTick', -0.6:0.3:0.6);
    % c= colorbar;
    % get(get(c, 'Children'))
    % c.Ticks = -0.6:0.1:0.6;
    if isempty(cmpath)
        % %{
        cmpath=['.\external\slanColor\自然配色\gallery\',num2str(5),'.jpg'];
        colorNum=256;
        CM = extract_colormap(cmpath, colorNum);
        colormap(flip(CM,1))
        % %}
        % colormap(plt.viridis)
        % colormap(flip(plt.viridis))
    end
else
    clim([0.05, 0.6]);
    % clm.Ticks = 0:0.3:0.6;
    % AxesH = axes('CLim', [0.05, 0.6]);
    % clm = colorbar('peer', AxesH, 'h', ...
    %            'XTick', 0:0.3:0.6);
    cbh = colorbar('h');
    get(get(cbh, 'Children'))
    if isempty(cmpath)
        % cmpath=['.\external\slanColor\自然配色\gallery\',num2str(4),'.jpg'];
        % colorNum=256;
        % CM = extract_colormap(cmpath, colorNum);
        % colormap(flip(CM,1))
        colormap(plt.viridis)
    end
end



% colormap(CM)

xlim([min(self.f, [], 'all'), max(self.f, [], 'all')])
ylim([min(self.f, [], 'all'), max(self.f, [], 'all')])

dtick = 5;
range = [0:dtick:max(self.f)];
hold on; plot_bsgrid(range)
title([])
colorbar off
shading interp
xticks(range)
yticks(range)
hA=gca;
hA.XAxis.FontSize=FontSize;
hA.YAxis.FontSize=FontSize;
xlabel('Frequency (Hz)')
ylabel('Frequency (Hz)')


end


% function plot_bc(self)
% 
% end

% 
% function plot_bchat(self)
% surfbc(self.fx,self.fy,abs(self.bchat));
% end



