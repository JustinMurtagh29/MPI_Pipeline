function gt = getFeaturesForEdges( p, gt )
% Get three different feature sets for given edges

segMeta = load([p.saveFolder 'segmentMeta.mat'], 'cubeIdx', 'segIds');
for i=1:length(gt)
    tic;
    uniqueCubes = unique(segMeta.cubeIdx(ismember(segMeta.segIds, gt(i).edges(:))));
    saveFolders = {p.local(uniqueCubes).saveFolder};
    featureFiles = cellfun(@(x)[x 'weights.mat'], saveFolders, 'uni', 0);
    [gt(i).oldFeatures, edgeIdx1] = getFeatures(gt(i).edges, saveFolders, featureFiles);
    featureFiles = cellfun(@(x)[x 'InterfaceFeatures.mat'], saveFolders, 'uni', 0);
    [gt(i).rawFeatures, edgeIdx2] = getFeatures(gt(i).edges, saveFolders, featureFiles);
    featureFiles = cellfun(@(x)[x 'InterfaceClassFeatures.mat'], saveFolders, 'uni', 0);
    [gt(i).classFeatures, edgeIdx3] = getFeatures(gt(i).edges, saveFolders, featureFiles);
    assert(all(edgeIdx1 == edgeIdx2) && all(edgeIdx2 == edgeIdx3));
    gt(i).edges = gt(i).edges(edgeIdx1,:);
    gt(i).labels = gt(i).labels(edgeIdx1);
    gt(i).prob = gt(i).prob(edgeIdx1);
    toc;
    display(num2str(i));
end

end

function [collectedFeatures, edgeIdx] = getFeatures(wantedEdges, saveFolders, featureFiles)

    collectedEdges = cell(size(saveFolders));
    collectedFeatures = cell(size(featureFiles));
    for i=1:length(saveFolders)
        % Load data and fix inconsitency in naming
        load([saveFolders{i} 'edges.mat']);
        load([saveFolders{i} 'borders.mat']);
        edges = edges(cat(1,borders(:).Area) > 10,:);
        load(featureFiles{i});
        if exist('weights', 'var')
            features = weights(cat(1,borders(:).Area) > 10,:);
            clear weights;
        elseif exist('X', 'var')
            features = X;
            clear X;
        end
        clear borders;
        % Keep only edges and features with one border as labels are ambigous otherwise
        [~, ~, idx] = unique(edges, 'rows');
        uniqueIdx = unique(idx);
        edgeHasMultipleBorder = uniqueIdx(histc(idx,uniqueIdx) > 1);
        edges(edgeHasMultipleBorder,:) = [];
        features(edgeHasMultipleBorder,:) = [];
        % Extract edges that were labelled in the dense tracing
        idx = ismember(edges, wantedEdges, 'rows');
        collectedEdges{i} = edges(idx,:);
        collectedFeatures{i} = features(idx,:);
    end
    collectedEdges = cat(1, collectedEdges{:});
    collectedFeatures = cat(1, collectedFeatures{:});
    [found, edgeIdx] = ismember(collectedEdges, wantedEdges, 'rows');
    assert(all(found));
    
end
