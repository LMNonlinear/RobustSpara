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
        imagePath = '.\external\slanColor\自然配色\gallery\4.jpg';  % Default to 256 colors, like turbo or hot in MATLAB
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
