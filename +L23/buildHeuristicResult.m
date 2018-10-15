% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/tmpscratch/amotta/l23/2018-10-09-mrnet-pipeline-run';
inDir = '/tmpscratch/amotta/l23/2018-10-13-vessel-and-nuclei-mask';
mag = [8, 8, 4];

vesselFile = fullfile(inDir, 'vessel_v1.mat');
nucleiFile = fullfile(inDir, 'nuclei_v1.mat');

info = Util.runInfo();
Util.showRunInfo(info);

%% Loading parameter
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

%% Run per-cube job
curTaskArgs = arrayfun( ...
    @(local) {local.bboxSmall}, ...
    param.local, 'UniformOutput', false);
curSharedArgs = {param.seg, vesselFile, nucleiFile, mag};

% Execute and retrieve results
job = Cluster.startJob( ...
    @jobFunction, curTaskArgs, ...
    'sharedInputs', curSharedArgs, ...
    'numOutputs', 3, 'cluster', {'memory', 12});
Cluster.waitForJob(job);
out = fetchOutputs(job);

%% Build output
result = table;
result.segIds = cat(1, out{:, 1});
result.vesselScore = cat(1, out{:, 2});
result.nucleiScore = cat(1, out{:, 3});

% NOTE(amotta): For backward compatibility
result.myelinScore = nan(size(vesselScore));

% Remove background and sort
result = result(result.segIds > 0, :);
result = sortrows(result, 'segIds');

result = table2struct(result, 'ToScalar', true);
result.info = info;

saveFile = fullfile(param.saveFolder, 'heuristicResult.mat');
Util.saveStruct(saveFile, result);

%% Main function
function [segIds, vesselScore, nucleiScore] = ...
        jobFunction(segParam, vesselFile, nucleiFile, mag, box)
    mag = reshape(mag, 1, []);
    magBox = box - box(:, 1) + 1;
    magBox = ceil(magBox ./ mag(:));
    magBox(:, 1) = magBox(:, 1) - mag(:) + 1;
    
    seg = loadSegDataGlobal(segParam, box);
   [segIds, ~, seg] = unique(seg(:));
    seg = reshape(seg, 1 + diff(box, 1, 2)');
    
    vesselScore = load(vesselFile);
    vesselScore = withMask(seg, mag, magBox, vesselScore.vesMask);
    
    nucleiScore = load(nucleiFile);
    nucleiScore = withMask(seg, mag, magBox, nucleiScore.nucMask);
end

function scores = withMask(seg, mag, magBox, mask)
    mask = mask( ...
        magBox(1, 1):magBox(1, 2), ...
        magBox(2, 1):magBox(2, 2), ...
        magBox(3, 1):magBox(3, 2));
    
    mask = double(mask);
    mask = imresize3(mask, size(mask) .* mag);
    scores = accumarray(seg(:), mask(:), [], @mean);
end
