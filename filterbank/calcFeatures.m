function calcFeatures(parameter, sub)
% We need to make sure that this function works correctly near borders of cube, right now due to no padding and 'same' convolution in filter3d
% we should have massive border effects?

% Why do this, empty initaliztion seems bad
weights = [];
segmentWeights = [];
weightNames = {};

% Create output directory if it does not exist
if ~exist(parameter.local(sub(1),sub(2),sub(3)).saveFolder, 'dir')
    mkdir(parameter.local(sub(1),sub(2),sub(3)).saveFolder);
end

% Load pixel list for borders and segments
load(parameter.local(sub(1),sub(2),sub(3)).borderFile);
load(parameter.local(sub(1),sub(2),sub(3)).segmentFile);
borders = borders'; % 1xN struct to Nx1 struct

% Define bounding box, add border so that 3D filter do not 
bbox = parameter.local(sub(1),sub(2),sub(3)).bboxSmall + [-10 10; -10 10; -10 10];

for l=1:length(parameter.feature.input)
    if strcmp(parameter.feature.input{l}, 'raw')
        imfeat = loadRawData(parameter.raw.root, parameter.raw.prefix, bbox, 1);
    elseif strcmp(parameter.feature.input{l}, 'aff')
        % Add Manuel: Here we need to differentiate between densly labeled regions and normal regions at least if we want to process LR beforehand (for training GP)
        if isfield(parameter.local(sub(1),sub(2),sub(3)), 'class')
            local = parameter.local(sub(1),sub(2),sub(3));
            imfeat = loadClassData(local.class.root, local.class.prefix, bbox);
        else
            imfeat = loadClassData(parameter.class.root, parameter.class.prefix, bbox);
        end
    end

    
    for m=1:size(parameter.filter,2)
        currentFilter = parameter.filter{m};
        for n=1:length(currentFilter{2})
            filteredImFeats = filter3d(currentFilter{1}, imfeat, currentFilter{2}(n), currentFilter{3});
            if isa(filteredImFeats, 'cell')
                for p=1:length(filteredImFeats)
                    filteredImFeats{p} = filteredImFeats{p}(11:end-10,11:end-10,11:end-10);
                    % Feature for edges
                    [weights_new, weightNames_new] = featureDesign(real(filteredImFeats{p}), borders);
                    weights = [weights weights_new];
                    % Feature for nodes
                    %[segmentWeights_new, weightNames_new] = featureDesign(real(filteredImFeats{p}), segments);
                    %segmentWeights = [segmentWeights segmentWeights_new]; 
                    % Names for features
                    weightNames_new = strcat({sprintf('%s %s filtersize=%d feature=%d ', currentFilter{1}, parameter.feature.input{l}, currentFilter{2}(n), p)}, weightNames_new);
                    weightNames = { weightNames{:}, weightNames_new{:} };
                end
            else
                filteredImFeats = filteredImFeats(11:end-10,11:end-10,11:end-10);
                % Feature for edges
                [weights_new, weightNames_new] = featureDesign(filteredImFeats, borders);
                weights = [weights weights_new];
                % Feature for nodes
                %[segmentWeights_new,weightNames_new] = featureDesign(filteredImFeats, segments);
                %segmentWeights = [segmentWeights segmentWeights_new];
                % Names for features
                weightNames_new = strcat({sprintf('%s %s filtersize=%d feature=%d ', currentFilter{1}, parameter.feature.input{l}, currentFilter{2}(n), 0)}, weightNames_new);
                weightNames = {weightNames{:}, weightNames_new{:}};                   
            end
        end
    end
end

% calculate shape features and add to weights
if isfield(parameter.local(sub(1),sub(2),sub(3)), 'class')
    siz = (parameter.local(sub(1),sub(2),sub(3)).bboxSmall(:,2) - parameter.local(sub(1),sub(2),sub(3)).bboxSmall(:,1) + 1)';
else
    siz = parameter.tileSize';
end

[weightsShape_borders, weightNames_new] = shapeFeatures(borders,siz);
%[weightsShape_segments, weightNames_new] = shapeFeatures(segments,siz);

weights = [weights weightsShape_borders];
%segmentWeights = [segmentWeights weightsShape_segments];
weightNames = {weightNames{:}, weightNames_new{:}};

save(parameter.local(sub(1),sub(2),sub(3)).weightFile, 'weights', 'weightNames');
%save(parameter.local(sub(1),sub(2),sub(3)).segmentWeightFile, 'segmentWeights', 'weightNames');

end

