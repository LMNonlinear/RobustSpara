% Author: Ying Wang, Min Li
% Create Time: 2025
% Copyright (c): 2020-2025 Ying Wang, yingwangrigel@gmail.com,
%                Min Li, minli.231314@gmail.com
% Joint China-Cuba LAB, UESTC, Hangzhou Dianzi University
% License: GNU General Public License v3.0 (see LICENSE file)

function CM = extract_colormap(imagePath, numColor)
    % Extracts a colormap from an image with interpolation to a specified number of colors.
    %
    % Parameters:
    % imagePath: string, the path to the image file.
    % colorNum: integer, the number of colors to extract (optional, default is 256).
    %
    % Returns:
    % CM: matrix, the interpolated colormap.
    
    % Set a default colorNum if not provided

    if nargin < 1 || isempty(imagePath)
        imagePath = '.\external\slanColor\鑷劧閰嶈壊\gallery\4.jpg';  % Default to 256 colors, like turbo or hot in MATLAB
    end
    if nargin < 2 || isempty(numColor)
        numColor = 256;  % Default to 256 colors, like turbo or hot in MATLAB
    end
    
    % Read the image from the specified path
    oriPic = imread(imagePath);
    
    % Convert the image to a list of RGB values
    % RGBList = double(reshape(oriPic, prod(size(oriPic, [1, 2])), 3));
    
    % Use rgb2ind to extract the colors
    [~, C] = rgb2ind(oriPic, 8);
    
    % Call the NTraveler function for color processing
    C = NTraveler(C);
    
    % Interpolate the colormap to 256 colors
    CM = interpColor(C, numColor);
    
    % Plot the colormap to visualize
    % figure;
    % colormap(CM);
    % colorbar;
    % title('Extracted Colormap');
end
