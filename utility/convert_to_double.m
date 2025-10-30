% Author: Ying Wang, Min Li
% Create Time: 2025
% Copyright (c): 2020-2025 Ying Wang, yingwangrigel@gmail.com,
%                Min Li, minli.231314@gmail.com
% Joint China-Cuba LAB, UESTC, Hangzhou Dianzi University
% License: GNU General Public License v3.0 (see LICENSE file)

function varargout = convert_to_double(varargin)
% Convert input data to double if it is not already
varargout = cell(1, nargin);
for i = 1:nargin

    if isa(varargin{i}, 'double')
        varargout{i} = varargin{i};

    else
        if isa(varargin{i}, 'gpuArray')
            varargout{i} = gather(varargin{i}); % Transfer data from GPU to MATLAB workspace
        else
            varargout{i} = varargin{i};
        end
        varargout{i} = double(varargout{i});
    end
end
end