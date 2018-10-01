function outputFile = makeSegmentationPreviewMovie(p, bbox, doClass, useClass)
% outputFile = makeSegmentationPreviewMovie(p, bbox, doClass)
%   Generates a segmentation preview movie (in AVI format)
%   with a given set of parameters.
%
%   p
%     Parameter structure
%
%   bbox
%     Bounding box around the region of interest. This can
%     be either a webKNOSSOS-style box (i.e., x_min, y_min,
%     z_min, x_max, y_max, z_mat) or MATLAB-style box (i.e.,
%     x_min, x_max; y_min, y_max; z_min, z_max).
%
%   doClass
%     This flag can be used to disable (and reuse) old cla-
%     ssification results for improved performance.
%     Default: true
%   
%   useClass
%   This flag can be used to create tempClass folder that
%   uses cnn output from p.class and writes to tempClass
%
%   outputFile
%     Path to generated segmentation preview movie

if ~exist('doClass', 'var')
    doClass = true;
end

if ~exist('useClass','var')
    useClass= false;
end

% also accept webKNOSSOS-style bouding boxes,
% i.e. (x_min, y_min, z_min, x_max, y_max, z_max)
if size(bbox, 1) == 1 || size(bbox, 2) == 1
    % if webKnossos style bounding box, transfer from width, height, depth to format used in pipeline
    bbox = Util.convertWebknossosToMatlabBbox(bbox);
end

% check size of bounding box
if ~Util.checkBoundingBox(bbox)
    error('Bounding box has invalid format');
end

% Classification (same function as bigFwdPass.m)
tempClass = struct;
tempClass.root = [p.tempFolder 'classForMovie/'];
tempClass.prefix = 'classForMovie';
tempClass.backend = 'wkwrap';

if doClass
    if isfield(tempClass, 'backend') && strcmp(tempClass.backend,'wkwrap')
        wkwInit('new', tempClass.root, 32, 32, 'single', 1);
    end

    classificationStep(p, tempClass, bbox);
else
    warning('> Skipping classification');
end

% create temp classification from p.class available output
if useClass
    if isfield(tempClass, 'backend') && strcmp(tempClass.backend,'wkwrap')
        wkwInit('new', tempClass.root, 32, 32, 'single', 1);
    end

    tempClassificationStep(p, tempClass, bbox);
    warning('> Using classification from p.class');
end

% Segmentation
tempSegFile = [p.tempFolder 'seg.mat'];
segmentationStep(p, tempClass, tempSegFile, bbox);

% Movie
outputFile = movieStep(p, tempSegFile, bbox);
end

function classificationStep(p, tempClass, bbox)
    functionH = @onlyFwdPass3DonKnossosFolder;
    inputCell{1} = {p.cnn.first, p.cnn.GPU, p.raw, tempClass, bbox, p.norm.func};
    jobName = 'classForMovie';

    if p.cnn.GPU
        job = startGPU(functionH, inputCell, jobName, 24);
    else
        job = startCPU(functionH, inputCell, jobName, 24);
    end

    Cluster.waitForJob(job);
end

function tempClassificationStep(p, tempClass, bbox)
data = loadClassData(p.class, bbox);
offset = bbox(:,1);
saveClassData(tempClass, offset, data);
end

function segmentationStep(p, tempClass, tempSegFile, bbox)
    functionH = @segmentForPipeline;
    inputCell{1} = {tempClass, bbox, p.seg.func, tempSegFile};
    jobName = 'segForMovie';

    job = startCPU(functionH, inputCell, jobName, 48);
    Cluster.waitForJob(job);
end

function outputFile = movieStep(p, tempSegFile, bbox)
    % Now low raw data
    raw = loadRawData(p.raw, bbox);
    raw = single(raw);

    % Normalize to [0 1] range as this is
    % what matlab expects of single images
    raw = raw ./ max(raw(:));

    % Load segmentation
    load(tempSegFile, 'seg');

    % Make movie
    outputFile = [p.tempFolder '_' num2str(p.seg.threshold) '_' datestr(clock, 30) '.avi']; % seg imposed on raw
    outputFile2 = [p.tempFolder '_' num2str(p.seg.threshold) '_' datestr(clock,30) '_RawAndSeg.avi']; % raw and seg side-by-side
    makeSegMovie(seg, raw, outputFile,outputFile2);
end
