function outputFile = makeSegmentationPreviewMovie(p, bbox, flag)
% Make movie to visualize segmentation given some parameter set p
if nargin <3
    flag = 1; % if this function is running in a loop this flag prevents it from calculating the classification multiple times
end

% also accept webKNOSSOS-style bouding boxes,
% i.e. (x_min, y_min, z_min, x_max, y_max, z_max)
if size(bbox, 1) == 1 || size(bbox, 2) == 1
    bbox = reshape(bbox, [3, 2]);
end

% First use same function as bigFwdPass.m to do classification
functionH = @onlyFwdPass3DonKnossosFolder;
tempClass.root = [p.tempFolder 'classForMovie/'];
tempClass.prefix = 'classForMovie';
inputCell{1} = {p.cnn.first, p.cnn.GPU, p.raw, tempClass, bbox, p.norm.func};

if flag == 1 % flag preventing the calculation classification for multiple movies
    if p.cnn.GPU
        job = startGPU(functionH, inputCell, 'classForMovie');
    else
        job = startCPU(functionH, inputCell, 'classForMovie');
    end
    Cluster.waitForJob(job);
else
    warning('omitting classification')
end

% Now do segmentation as in miniSegmentation.m
tempSegFile = [p.tempFolder 'seg.mat'];
inputCell{1} = {tempClass.root, tempClass.prefix, bbox, p.seg.func, tempSegFile};
functionH = @segmentForPipeline;
job = startCPU(functionH, inputCell, 'segForMovie');
Cluster.waitForJob(job);

% Now low raw data
raw = loadRawData(p.raw.root, p.raw.prefix, bbox, false);

% Normalize to [0 1] range as this is what matlab expects of single images
raw = raw ./ max(raw(:));

% Load segmentation
load(tempSegFile);

% Make movie
outputFile = [p.tempFolder datestr(clock, 30) '.avi'];
makeSegMovie(seg, raw, outputFile);

end

