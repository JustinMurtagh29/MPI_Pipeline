function job = calculateFeaturesForTrainingData( dataFolder, ...
    saveFolder, featureMap, ignoreBorder, cluster, discardUndecided )
%CALCULATEFEATURESFORTRAININGDATA Calculate the specified feature map for
%all training cube in dataFolder.
% INPUT dataFolder: string
%           Path to training data folder.
%       saveFolder: string
%           Path to folder where features are saved.
%       featureMap: Interface.FeautreMap object
%       ignoreBorder: (Optional) logical
%           Flat to set the feature map border to zero during calculation
%           (compatibility with old synapse detection data).
%           (Default: false)
%       cluster: (Optional) parallel.cluster object
%           Cluster object used for calculation.
%           (Default: getCluster('cpu') if this function exists)
%       discardUndecided: (Optional) logical
%           Flag indicating that interfaces annotated with undecided are to
%           be discarded.
%           (Default: false)
% OUTPUT job: job object
%           Job object for feature calculation on workers.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if ~exist('ignoreBorder','var') || isempty(ignoreBorder)
    ignoreBorder = false;
end
if ~exist('discardUndecided','var') || isempty(discardUndecided)
    discardUndecided = false;
end
subvols = find(ismember([40, 80, 160], featureMap.subvolsSize));

s = what(dataFolder);
inputCell = cell(length(s.mat),1);
for i = 1:length(s.mat)
    inputCell{i} = {[SynEM.Util.addFilesep(dataFolder), s.mat{i}], ...
        SynEM.Util.addFilesep(saveFolder), featureMap, ignoreBorder, ...
        subvols, discardUndecided};
end


try
    if ~exist('cluster','var') || isempty(cluster)
        cluster = getCluster('cpu');
    end
    job = SynEM.Util.startJob(cluster,@jobWrapper,inputCell,0);
catch
    warning(['No cluster object found. ' ...
        ' Calculation is done sequentially.']);
    cellfun(@(x)jobWrapper(x{:}), inputCell);
    job = [];
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function jobWrapper(file, saveFolder, featureMap, useZeroBorder, ...
    subvols, discardUndecided)

Util.log('Processing file %s', file);

% featureMap.featTexture{2}.options.gpuDev = true;

%load data
m = load(file);
data = m.data;
raw = single(data.raw);
interfaceLabels = m.interfaceLabels;
undecidedList = m.undecidedList;

% discard undecided interfaces
if discardUndecided
    Util.log('Discarding undecided interfaces.');
    data.interfaceSurfaceList(undecidedList) = [];
    data.subsegmentsList = cellfun(@(x)x(~undecidedList,:), ...
        data.subsegmentsList, 'UniformOutput', false);
    interfaceLabels(undecidedList) = [];
end

% discard interfaces at border
if isfield(m ,'borderInt')
    Util.log('Discarding interfaces at border.');
    isAtBorder = m.borderInt.isAtBorder;
    if discardUndecided
        isAtBorder(undecidedList) = [];
    end
    data.interfaceSurfaceList(isAtBorder) = [];
    data.subsegmentsList = cellfun(@(x)x(~isAtBorder,:), ...
        data.subsegmentsList, 'UniformOutput', false);
    interfaceLabels(isAtBorder) = [];
end

%apply area threshold and select subvols (required for labels)
areaT = cellfun(@length, data.interfaceSurfaceList) > featureMap.areaT;
interfaces.surface = data.interfaceSurfaceList(areaT);
interfaces.subseg = cellfun(@(x)x(areaT,:),data.subsegmentsList, ...
    'UniformOutput',false);
interfaces.subseg = interfaces.subseg(subvols);
interfaceLabels = interfaceLabels(areaT);

%calculate outputs
X = featureMap.calculate(interfaces, raw, useZeroBorder);
y = [interfaceLabels == 1; interfaceLabels == 2];

%save results
if ~exist(saveFolder,'dir')
    mkdir(saveFolder)
end
[~,name] = fileparts(file);
Util.log('Saving result to %s.', [saveFolder, 'Features_', name, '.mat']);
m = matfile([saveFolder, 'Features_', name, '.mat'],'Writable',true);
m.X = X;
m.classLabels = y;
m.featureMap = featureMap;

end
