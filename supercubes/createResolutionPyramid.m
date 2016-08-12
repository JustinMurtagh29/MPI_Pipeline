function createResolutionPyramid( root, prefix, bbox, outputDirectory, subsamplingSegmentation )
%Input:
%   If using just one input argument, and normal wK hierachy with JSONs already present, rest inferred
%   root = Directory of KNOSSOS Hierachy mag1
%   prefix = File prefix of KNOSSOS Hierachy mag1
%   bbox = Bounding Box of KNOSSOS Hierachy mag1 
%       (minimal coordinates should align with start of Knossos cube (e.g. [1 1 1]))
%       (Not sure whether it will work if lower corner in bbox is not [1 1 1], NOT tested)
%   outputDirectory = where to write higher resolutions (subfolder named as
%       the magnification will be created inside)
%   subsamplingSegmentation = if you try to subsample a segmentation (instead of raw data)
%       set this to true, will use @mode instead of @median for subsampling

if nargin < 5
    subsamplingSegmentation = false;
end

if nargin < 4
    % If no output folder is provided, assume one level up from root
    fsep = filesep;
    outputDirectory = strrep(root, [fsep '1' fsep], fsep);
    clear fsep;
end

if nargin < 3
    % Get bounding box from section.json
    temp = readJson([outputDirectory 'section.json']);
    bbox = double(cell2mat(cat(2, temp.bbox{:}))');
    % Make sure it does not start at [0 0 0] cause that would make x-0001 Knossos hierachy folder to be written
    % wK vs. matlab coordinate system offset seems to create some confusion (as it should)
    bbox = bsxfun(@max, [1; 1; 1], bbox);
    clear temp;
end

if nargin < 2
    % Get prefix from filename in root Knossos Hierachy
    temp = dir([root 'x0000/y0000/z0000/*.raw']);
    [startIdx, endIdx] = regexp(temp(1).name, '_mag.{1,3}_');
    prefix = temp(1).name(1:endIdx-1);
    clear temp startIdx endIdx;
end

% Create output folder if it does not exist
if ~exist(outputDirectory, 'dir');
    mkdir(outputDirectory);
end

% Set according to memory limits, currently optimized for 12 GB RAM, segmentation will need 48 GB currently
% Will be paralellized on cubes of this size
cubeSize = [1024 1024 1024];
% Write these magnifications, mags grouped in inside Brackets will be calculated by only reading once
magsToWrite = {[2 4 8] [16 32 64] [128 256 512]};
% set function for downsampling based on whether subsampling raw data or segmentation
if subsamplingSegmentation
    downsamplingFunction = @mode;
else
    downsamplingFunction = @median;
end

% Do the work, submitted to scheduler
bboxTemp = bbox;
functionH = @writeSupercubes;
for i=1:length(magsToWrite)
    X = unique([bboxTemp(1,1):cubeSize(1):bboxTemp(1,2) bboxTemp(1,2)]);
    Y = unique([bboxTemp(2,1):cubeSize(2):bboxTemp(2,2) bboxTemp(2,2)]);
    Z = unique([bboxTemp(3,1):cubeSize(3):bboxTemp(3,2) bboxTemp(3,2)]);
    idx = 1;
    for x=1:length(X)-1
        for y=1:length(Y)-1
            for z=1:length(Z)-1
                thisBBox = [X(x) X(x+1)-1; Y(y) Y(y+1)-1; Z(z) Z(z+1)-1];
                inputCell{idx} = {root, prefix, thisBBox, magsToWrite{i}, outputDirectory, subsamplingSegmentation};
                idx = idx + 1;
            end
        end
    end
    job = startCPU(functionH, inputCell, 'supercubes');
    wait(job, 'finished');
    % Update parameter, read last resoultion written with right bbox
    bboxTemp = ceil((bbox - 1)./ magsToWrite{i}(end)) + 1;
    root = fullfile(outputDirectory, num2str(magsToWrite{i}(end)));
    [startIdx, endIdx] = regexp(prefix, '_mag.{1,3}');
    prefix = [prefix(1:startIdx-1) '_mag' num2str(magsToWrite{i}(end))];
    clear inputCell; 
end

end

