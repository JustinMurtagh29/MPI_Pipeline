function calcFeatures(parameter, sub)
% We need to make sure that this function works correctly near borders of cube, right now due to no padding and 'same' convolution in filter3d
% we should have massive border effects?

weights = [];

if ~exist(parameter.local(sub(1),sub(2),sub(3)).saveFolder, 'dir')
    mkdir(parameter.local(sub(1),sub(2),sub(3)).saveFolder);
end

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
end

load(parameter.local(sub(1),sub(2),sub(3)).borderFile);
for m=1:size(parameter.filter,2)
    currentFilter = parameter.filter{m};
    for n=1:length(currentFilter{2})
        filteredImFeats = filter3d(currentFilter{1}, imfeat, currentFilter{2}(n), currentFilter{3});
        if isa(filteredImFeats, 'cell')
            for p=1:length(filteredImFeats)
                filteredImFeats{p} = filteredImFeats{p}(11:end-10,11:end-10,11:end-10);
                weights_new = featureDesign(real(filteredImFeats{p}), borders);
                weights = [weights weights_new];
                disp(sprintf('%s size %d feature %d x%d', currentFilter{1}, currentFilter{2}(n), p, size(weights_new, 2)));
            end
        else
            filteredImFeats = filteredImFeats(11:end-10,11:end-10,11:end-10);
            weights_new = featureDesign(filteredImFeats, borders);
            weights = [weights weights_new];
            disp(sprintf('%s size %d feature %d x%d', currentFilter{1}, currentFilter{2}(n), 0, size(weights_new, 2)));
        end
    end
end

% calculate shape features and add to weights
if isfield(parameter.local(sub(1),sub(2),sub(3)), 'class')
    siz = (parameter.local(sub(1),sub(2),sub(3)).bboxSmall(:,2) - parameter.local(sub(1),sub(2),sub(3)).bboxSmall(:,1) + 1)';
else
    siz = parameter.tileSize';
end

weightsShape_borders = shapeFeatures(borders,siz);

weights = [weights weightsShape_borders];

save(parameter.local(sub(1),sub(2),sub(3)).weightFile, 'weights');

end

