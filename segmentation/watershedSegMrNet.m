function seg = watershedSegMrNet(class, varargin)
    % seg = watershedSegMrNet(class)
    %   Generates a watershed-based over-segmentation using log-distance
    %   transformed results from Benedikt's MrNet.
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    opt = struct;
    opt.classThresh = 0;
    opt.minDepth = 0.2;
    opt.voxelSize = [1, 1, 1];
    opt = Util.modifyStruct(opt, varargin{:});
    
    reSample = round(opt.voxelSize ./ min(opt.voxelSize));
    assert(all(reSample(1:2) == 1));
    
    %% Distance transform
    dist = class < opt.classThresh;
    dist = repelem(dist, 1, 1, reSample(3));
    dist = -log(1 + bwdist(dist));
    dist = dist(:, :, 1:reSample(3):end);
    
    %% Define seed locations
    dist = imhmin(dist, opt.minDepth, 26);
    minima = imregionalmin(dist, 26);
    dist = imimposemin(dist, minima, 26);
    
    %% Perform watershed
    seg = watershed(dist, 26);
end

