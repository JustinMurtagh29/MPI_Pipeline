function X = calculate(obj, interfaces, raw, ignoreBorder)
%CALCULATE The main entry point for the feature calculation for interface
% classification. This function iterates over all features in the feature
% map and calculates the selected features.
% INPUT interfaces: Struct containing the fields
%           'surface': [Nx1] cell array. Each cell contains the linear
%               indices of an interface/border between two segments.
%           'subseg': [1xM] cell array where M = obj.numSubvolumes. Each
%               cell contains a [Nx2] cell array where N is the same as for
%               the 'surface' field. Each cell in this array contains the
%               linear indices of a subsegment of the corresponding
%               subvolume for the interface.
%       raw: 3d numerical array containing the raw data. Will be cast to
%           single during calculation before normalization. The size of raw
%           should be the size of the cube to which the indices in
%           interfaces refer plus the border of the feature map, i.e. each
%           dimension is cropped by border./2 on each side of raw.
%           (see also obj.fRawNorm)
%       ignoreBorder: (Optional) Logical specifying whether feature maps
%           from texture feature calculation are cropped by the interface
%           border.
%           (Default: false).
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

%compatibility with old data
if ~exist('ignoreBorder','var') || isempty(ignoreBorder)
    ignoreBorder = false;
end

%check for number of subvolumes
if obj.numSubvolumes ~= length(interfaces.subseg)
    error(['Number of subsgments in interfaces does not equal the ' ...
        'number of subsegments specified in the feature map.']);
end

%texture features
if isempty(obj.selectedFeat)
    obj.setSelectedFeat(true(obj.numFeatures,1));
end

%apply area threshold and prepare pooling volumes
indices = cellfun(@(x)length(x) > obj.areaT,interfaces.surface);
interfaces.surface = interfaces.surface(indices);
if obj.numSubvolumes > 0
    interfaces.subseg = cellfun(@(x)x(indices,:), ...
        interfaces.subseg,'UniformOutput',false);
    subsegmentsList = interfaces.subseg;
    subsegmentsList = cat(2,subsegmentsList{1:length(subsegmentsList)});
    vols = cat(2,interfaces.surface, subsegmentsList);
else
    vols = interfaces.surface;
end
fNorm = str2func(obj.fRawNorm);
raw = fNorm(raw);

%preallocate output
switch obj.mode
    case 'direction'
        X = zeros(2*length(interfaces.surface), ...
            sum(obj.numFeaturesSelected),'like', raw);
    case 'single'
        X = zeros(length(interfaces.surface), ...
            sum(obj.numFeaturesSelected),'like', raw);
end


%calculate texture features
count = 1;
for i = 1:length(obj.featTexture)
    numFeat = sum(obj.selectedFeat{i}(:));
    if numFeat > 0
        fprintf(['[%s] SynEM.FeatureMap.calculate - Calculating ' ...
            'texture feature %d: %s.\n'], datestr(now), i, ...
            obj.featTexture{i}.name);
        X(1:length(interfaces.surface),count: count + numFeat - 1) = ...
            obj.calculateTextureFeature( vols, raw, i, ignoreBorder);
        count = count + numFeat;
    end
end

%reference cube size
if ignoreBorder
    cubeSize = size(raw);
else
    cubeSize = size(raw) - obj.border;
end

%transform to 3d coordinates once for all voxel features
if ~isempty(obj.featShape) %only do it if necessary
    if obj.numSubvolumes > 0
        coordVols = cellfun(@(x)obj.ind2sub(cubeSize,x), ...
            vols(:,[1 end-1 end]),'UniformOutput',false);
    else
        coordVols = cellfun(@(x)obj.ind2sub(cubeSize,x), ...
            vols(:,1),'UniformOutput',false);
    end
end

%calculate shape features.
for i = 1:length(obj.featShape)
    numFeat = sum(obj.selectedFeat{i + length(obj.featTexture)});
    if numFeat > 0
        fprintf(['[%s] SynEM.FeatureMap.calculate - Calculating shape ' ...
            'feature %i: %s.\n'], datestr(now), i, obj.featShape{i}.name);
        X(1:length(interfaces.surface),count:count + numFeat - 1) = ...
            obj.calculateShapeFeatures(i, coordVols);
        count = count + numFeat;
    end
end

%direction inversion of interfaces
if strcmp(obj.mode, 'direction')
    X((length(interfaces.surface) + 1): end,:) = obj.invertDirection( ...
        X(1:length(interfaces.surface),:));
end

if any(isnan(X(:)))
    warning('Feature matrix contains NaNs.')
end
end
