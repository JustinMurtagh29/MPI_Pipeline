function show(param, cubeIdx)
    % show(param, cubeIdx)
    %   Show a movie with the raw data and overlaid colours
    %   for a subset of the borders in the specified cube.
    %   This function is intended for debugging purposes.
    %
    % param
    %   Parameter structure produced by setParameterSettings
    %
    % cubeIdx
    %   Linear index of the cube of interest
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    % config
    pickCount = 100;

    cubeParam = param.local(cubeIdx);
    cubeDir = cubeParam.saveFolder;
    
    % load extended borders
    data = load([cubeDir, 'bordersExt.mat']);
    
    % extract data
    box = data.box;
    borders = data.borders;
    
    disp('Loading raw data...');
    raw = loadRawData(param.raw, box);
    
    % picking a subset of edges
    edgeCount = size(borders, 1);
    edgeIds = randperm( ...
        edgeCount, min(edgeCount, pickCount));
    
    % highlight borders and show result
    colored = highlightBorders( ...
        raw, 0.5, borders(edgeIds, :));
    implay(colored);
end

function data = gray2rgb(data)
    dataLen = size(data, 3);
    data = arrayfun(@(idx) ...
        repmat(data(:, :, idx), [1, 1, 1, 3]), ...
        1:dataLen, 'UniformOutput', false);
    data = cat(3, data{:});
end

function data = highlightBorders(data, alpha, borders)
    % config
    colors = [ ...
          0,   0,   0;  % black
        255, 255,   0;  % yellow
        255,   0,   0]; % red
    colors = uint8(colors);
    
    % NOTE
    %   MATLAB represents a color movie as a KxLx3xN matrix
    %   where KxL repesent a XY slice, the third dimension
    %   is the color channels (RGB) and the fourth dimension
    %   make up the different frames.
    %
    %   For processin and coloring in borders, it's easier to
    %   work with a KxLxNx3 matrix. Hence the permutation of
    %   dimensions at the end.
    
    % mix in colors
    edgeCount = size(borders, 1);
    
    % convert data to RGB
    dataCount = numel(data);
    data = gray2rgb(data);
    
    for curIdx = 1:edgeCount
        for curSubIdx = 1:3
            % get voxel indices
            curVoxelIds = borders{curIdx, curSubIdx};
            
            % iterate over color channels
            for curRGB = 1:3
                curRgbVoxelIds =  ...
                    curVoxelIds + (curRGB - 1) * dataCount;
                data(curRgbVoxelIds) =  ...
                    (1 - alpha) * data(curRgbVoxelIds) ...
                  + alpha * colors(curSubIdx, curRGB);
            end
        end
    end
    
    % fix dimensions
    data = permute(data, [1, 2, 4, 3]);
end