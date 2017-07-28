function job = cnnPred( cnet, dataFolder, saveFolder, p, options, cluster )
%CNNPRED Calculate CNN predictions for the training data cubes.
% INPUT cnet: Codat.CNN.cnn object or string
%           CNN for prediction. If cnet is a string the it is interpreted
%           as the path to a matfile containing the cnet in a variable
%           called cnet.
%       dataFolder: string
%           Path to training data folder.
%       saveFolder: string
%           Path to folder where features are saved.
%       p: struct
%           Segmentation parameter struct to load the raw data.
%       options: struct
%           Options for prediction (see Codat.CNN.Misc.predictCube).
%       cluster: (Optional) parallel.cluster object
%           Cluster object for parallel calculation.
%           (Default: Calculation is done sequentially)
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

%create output folder
saveFolder = Util.addFilesep(saveFolder);
if ~exist('saveFolder', 'dir')
    mkdir(saveFolder)
end

%store input parameters
paramsFile = fullfile(saveFolder, 'InputParams.mat');
Util.log('Saving parameters to %s.', paramsFile);
git_hash = Util.getGitHash();
Util.save(paramsFile, cnet, git_hash, dataFolder, saveFolder, options);

if ischar(cnet)
    [~,cname] = fileparts(cnet);
    copyfile(cnet, fullfile(saveFolder, [cname '.mat']));
    m = load(cnet);
    cnet = m.cnet;
end

if ~exist('options', 'var')
    options = struct;
    options.target_size = [128, 128, 128];
    options.display = 10;
end

%prepare worker data
s = what(dataFolder);
inputCell = cell(length(s.mat),1);
for i = 1:length(s.mat)
    inputCell{i} = {cnet, fullfile(dataFolder, s.mat{i}), saveFolder, ...
        options, p};
end

if exist('cluster', 'var') && ~isempty(cluster)
    Util.log('Submitting jobs to cluster.');
    job = Cluster.startJob(@jobWrapper, inputCell, 'cluster', cluster);
else
    Util.log('No cluster input found. Calculation is done sequentially.');
    for i = 1:length(inputCell)
        jobWrapper(inputCell{i});
    end
    job = [];
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function jobWrapper(cnet, dataFile, saveFolder, options, p)
Util.log('Prediction for file %s.', dataFile);
m = load(dataFile);
bbox = m.metadata.bboxBig;
raw = Seg.IO.loadRaw(p, bbox, cnet.border./2);
raw = (single(raw) - 122)./22;

pred = Codat.CNN.Misc.predictCube(cnet, raw, options);

[~,name] = fileparts(dataFile);
saveName = [saveFolder, 'Features_', name, '.mat'];
Util.log('Saving result to %s.', saveName);
Util.save(saveName, pred, bbox);
end

