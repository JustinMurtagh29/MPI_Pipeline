function gt = getFeaturesForEdges( p, gt, rawFeatFile, classFeatFile )
% Get three different feature sets for given edges

if ~exist('rawFeatFile', 'var') || isempty(rawFeatFile)
    rawFeatFile = 'InterfaceRawFeatures.mat';
end

if ~exist('classFeatFile', 'var') || isempty(classFeatFile)
    classFeatFile = 'InterfaceClassFeatures.mat';
end

segMeta = load([p.saveFolder 'segmentMeta.mat'], 'cubeIdx', 'segIds');
for i=1:length(gt)
    tic;
    uniqueCubes = unique(segMeta.cubeIdx(ismember(segMeta.segIds, gt(i).edges(:))));
    saveFolders = {p.local(uniqueCubes).saveFolder};
    bboxes = {p.local(uniqueCubes).bboxSmall};
    featureFiles = cellfun(@(x)[x rawFeatFile], saveFolders, 'uni', 0);
    [gt(i).rawFeatures, edgeIdx1, gt(i).borderCoM] = getFeatures(gt(i).edges, saveFolders, featureFiles, bboxes);
    featureFiles = cellfun(@(x)[x classFeatFile], saveFolders, 'uni', 0);
    [gt(i).classFeatures, edgeIdx2, test1] = getFeatures(gt(i).edges, saveFolders, featureFiles, bboxes);
    assert(all(edgeIdx1 == edgeIdx2));
    assert(all(gt(i).borderCoM(:) == test1(:)));
    gt(i).edges = gt(i).edges(edgeIdx1,:);
    gt(i).labels = gt(i).labels(edgeIdx1);
    gt(i).prob = gt(i).prob(edgeIdx1);
    toc;
    display(num2str(i));
end

end

function [collectedFeatures, edgeIdx, collectedBorderCoM] = getFeatures(wantedEdges, saveFolders, featureFiles, bboxes)

    collectedEdges = cell(size(saveFolders));
    collectedFeatures = cell(size(featureFiles));
    for i=1:length(saveFolders)
        % Load data
        load([saveFolders{i} 'edges.mat']);
        load([saveFolders{i} 'borders.mat']);
        edges = edges(cat(1,borders(:).Area) > 10,:);
        borderCoM = round(bsxfun(@plus, cat(1,borders(:).Centroid), bboxes{i}(:,1)' - 1));
        borderCoM = borderCoM(cat(1,borders(:).Area) > 10,:);
        clear borders;
        load(featureFiles{i});
        % Not needed anymore due to focused annotation
%         % Keep only edges and features with one border as labels are ambigous otherwise
%         [~, ~, idx] = unique(edges, 'rows');
%         uniqueIdx = unique(idx);
%         edgeHasMultipleBorder = uniqueIdx(histc(idx,uniqueIdx) > 1);
%         edges(ismember(idx, edgeHasMultipleBorder),:) = [];
%         features(ismember(idx, edgeHasMultipleBorder),:) = [];
        % Extract edges that were labelled in the dense tracing
        idx = ismember(edges, wantedEdges, 'rows');
        collectedEdges{i} = edges(idx,:);
        collectedFeatures{i} = features(idx,:);
        collectedBorderCoM{i} = borderCoM(idx,:);
    end
    collectedEdges = cat(1, collectedEdges{:});
    collectedFeatures = cat(1, collectedFeatures{:});
    collectedBorderCoM = cat(1, collectedBorderCoM{:});
    [found, edgeIdx] = ismember(collectedEdges, wantedEdges, 'rows');
    assert(all(found));
    
end
