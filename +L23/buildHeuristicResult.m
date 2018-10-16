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
curSharedArgs = {param, vesselFile, nucleiFile, mag};

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
result.myelinScore = nan(size(result.vesselScore));

% Remove background and sort
result = result(result.segIds > 0, :);
result = sortrows(result, 'segIds');

result = table2struct(result, 'ToScalar', true);
result.info = info;

saveFile = fullfile(param.saveFolder, 'heuristicResult.mat');
Util.saveStruct(saveFile, result);

%% Main function
function [segIds, vesselScore, nucleiScore] = ...
        jobFunction(param, vesselFile, nucleiFile, mag, box)
    mag = reshape(mag, 1, []);
    
    seg = loadSegDataGlobal(param.seg, box);
   [segIds, ~, seg] = unique(seg(:));
    seg = reshape(seg, 1 + diff(box, 1, 2)');
    
    vesselScore = load(vesselFile);
    vesselScore = withMask(param, mag, box, seg, vesselScore.vesMask);
    
    nucleiScore = load(nucleiFile);
    nucleiScore = withMask(param, mag, box, seg, nucleiScore.nucMask);
end

function scores = withMask(param, mag, box, seg, mask)
    magBox = box - param.bbox(:, 1) + 1;
    magBox = ceil(magBox ./ mag(:));
    
    maskOff = (magBox - [1, 0]) .* mag(:) - [0, 1];
    maskOff = maskOff + param.bbox(:, 1) - box;
    
    mask = mask( ...
        magBox(1, 1):magBox(1, 2), ...
        magBox(2, 1):magBox(2, 2), ...
        magBox(3, 1):magBox(3, 2));
    
    mask = double(mask);
    mask = imresize3(mask, size(mask) .* mag, 'linear');
    
    mask = mask( ...
        (1 + maskOff(1, 1)):(end - maskOff(1, 2)), ...
        (1 + maskOff(2, 1)):(end - maskOff(2, 2)), ...
        (1 + maskOff(3, 1)):(end - maskOff(3, 2)));
    
    scores = accumarray(seg(:), mask(:), [], @mean);
end
