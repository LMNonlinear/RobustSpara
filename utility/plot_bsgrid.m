function plot_bsgrid(range)
% range: The maximum value for the x:axis (e.g., max(self.f))

% Sample data
% x = range;
% y = x;  % Example line with :1 slope, modify as needed
hold on;
% c=[0.3,0.3,1];
c = [1, 0, 0];
range_f1_f2 = [range(:); range(:) + max(range(:)) + range(2)-range(1)]';
%
% % Add vertical lines at every integer point from 0 to range
xline(range, 'Color', c, 'LineStyle', ':', 'LineWidth', 0.2);
yline(range, 'Color', c, 'LineStyle', ':', 'LineWidth', 0.2);
% % Add lines with :1 slope
% for i = 1:range
%     xline(i, 'b:');  % Modify color/style as needed
%     yline(i, 'b:');  % Modify color/style as needed
% end

% % Optional: Add labels to the axes
% xlabel('X:axis');
% ylabel('Y:axis');
% Plot the data
% z=mean(zlim);
z = (zlim);
z = z(2);

for i = range_f1_f2
    h = plot3([0, i], [i, 0], [z, z], 'Color', c, 'LineStyle', ':', 'LineWidth', 0.2); %,'Layer', 'top' zorder
    % h.Visible="on";
    % g.Children
    uistack(h, 'top');
end

end
