function gt = getContinuityLabelsFromNml(p, pT, fixOffset)
% Extract labels for interface continuity classification from a set of nml
% files

for i=1:length(pT.local)
    display(['----- Training region ' num2str(i) ' -----']);
    for j=1:length(pT.local(i).trainFile)
        display(['----- Training file ' num2str(j) ' -----']);
        % Load & process dense merger mode skeleton files
        skel = skeleton(pT.local(i).trainFile{j});
        % Extract node positions and corresponding segment ids
        nodes = cellfun(@(x)x(:,1:3), skel.nodes, 'uni', 0);   
        nodesPerTree = cellfun(@(x)size(x,1), nodes);
        nodes = cat(1,nodes{:});
        
        if fixOffset
            % Temporary fix: Empirical fix to not having to reuse ex145-conversion
            nodes = bsxfun(@minus, nodes, [2 2 1]);
        end
        
        segIdsOfGT = Seg.Global.getSegIds(p, nodes);
        segIdsOfGT = mat2cell(segIdsOfGT, nodesPerTree);
        % Remove nodes placed in background (0)
        nrNodesInBackground = sum(cellfun(@(x)sum(x==0), segIdsOfGT));
        display(['Removing ' num2str(nrNodesInBackground) ' of ' ...
            num2str(sum(nodesPerTree)) ' nodes due to background placement!']);
        segIdsOfGT = cellfun(@(x)x(x~=0), segIdsOfGT, 'uni', 0);
        % Keep only unique hits per tree
        segIdsOfGT = cellfun(@(x)unique(x), segIdsOfGT, 'uni', 0);

        % Load segmentation data to determine all segments in bounding box
        seg = loadSegDataGlobal(p.seg, pT.local(i).bboxSmall);
        segMeta = load([p.saveFolder 'segmentMeta.mat'], 'cubeIdx', 'segIds', 'voxelCount');
        segIdsInBBox = unique(seg(seg~=0));
        clear seg;

        % Remove merged segments from list
        allSegmentIdsInGT = unique(cat(1,segIdsOfGT{:}));
        segmentsCounts = histc(cat(1,segIdsOfGT{:}), allSegmentIdsInGT);
        excludedSegments = allSegmentIdsInGT(segmentsCounts > 1);
        display(['Removing ' num2str(numel(excludedSegments)) ' of ' ...
            num2str(sum(nodesPerTree)) ' segments due to being covered by more than one tree!']);
        segIdsOfGT = cellfun(@(x)x(~ismember(x, excludedSegments)), segIdsOfGT, 'uni', 0);
        clear segmentIdsInGT;

        % Remove segments not in bounding box
        segIdsOfGT = cellfun(@(x)x(ismember(x,segIdsInBBox)), segIdsOfGT, 'uni', 0);

        % Get all edges between segments in bounding box
        uniqueCubes = unique(segMeta.cubeIdx(ismember(segMeta.segIds, segIdsInBBox)));
        [edges, prob] = Seg.Global.getGlobalEdges(p, uniqueCubes);
        edgeIdxInBbox = all(ismember(edges, segIdsInBBox),2);
        edgesInBbox = edges(edgeIdxInBbox,:);
        probInBbox = prob(edgeIdxInBbox,:);
        clear edgeIdxInBbox edges uniqueCubes;

        % Determine labels
        labels = determineLabelsFromEqClass(segIdsOfGT, edgesInBbox);

        % Calculate and output some statistics
        statisticsOfDenseTracing(segMeta, segIdsInBBox, segIdsOfGT)

        % Some statistics as well
        statisticsOfLabels(labels);

        % Transfer to structure with clear naming
        gt(i,j).edges = edgesInBbox;
        gt(i,j).prob = probInBbox;
        gt(i,j).labels = labels;
    end
end

end

function statisticsOfDenseTracing(segMeta, segIdsInBBox,segIdsOfGT)

    % Some statistics output
    segmentIdsInGT = unique(cat(1,segIdsOfGT{:}));
    segmentsCounts = histc(cat(1,segIdsOfGT{:}), segmentIdsInGT);
    nrSegmentsInBoundingBox = numel(segIdsInBBox);
    nrSegmentsInSkeletons = numel(segmentIdsInGT);
    volSegmentsInBoundingBox = sum(segMeta.voxelCount(segIdsInBBox));
    volSegmentsInSkeletons = sum(segMeta.voxelCount(segmentIdsInGT));

    % Print to console for now
    display(['Segments in bounding box: ' ...
        num2str(nrSegmentsInBoundingBox)]);
    display(['Segments collected by dense tracing: ' ...
        num2str(nrSegmentsInSkeletons)]);
    display(['Segments that are part of more than one tree: ' ...
        num2str(sum(segmentsCounts > 1))]);
    display(['Fraction collected segments: ' ...
        num2str(nrSegmentsInSkeletons./nrSegmentsInBoundingBox, '%3.2f')]);
    display(['Fraction collected volume: ' ...
        num2str(volSegmentsInSkeletons./volSegmentsInBoundingBox, '%3.2f')]);
    
end

function statisticsOfLabels(labels)

    % Print to console for now
    display(['Positive labels: ' ...
        num2str(sum(labels == 1)) ', ' num2str(sum(labels == 1)./numel(labels)*100) ' %']);
    display(['Negative labels: ' ...
        num2str(sum(labels == -1)) ', ' num2str(sum(labels == -1)./numel(labels)*100) ' %']);
    display(['Unlabeled edges: ' ...
        num2str(sum(labels == 0)) ', ' num2str(sum(labels == 0)./numel(labels)*100) ' %']);
    
end

function labels = determineLabelsFromEqClass(segIdsOfGT, edgesInBbox)
    labels = zeros(size(edgesInBbox,1),1, 'int8');
    temp = segIdsOfGT(cellfun(@numel, segIdsOfGT) > 1);
    edgesPos = cellfun(@(x)combnk(x,2), temp, 'uni', 0);
    edgesPos = cat(1, edgesPos{:});
    edgesNeg = cell(size(segIdsOfGT));
    for i=1:length(segIdsOfGT)
        otherSegments = cat(1, segIdsOfGT{setdiff(1:length(segIdsOfGT), i)});
        theseSegments = segIdsOfGT{i};
        edgesNeg{i} = cat(2, repmat(otherSegments, size(theseSegments)), ...
            repmat(theseSegments, size(otherSegments)));
    end
    edgesNeg = cat(1, edgesNeg{:});
    edgesNeg = sort(edgesNeg, 2);
    assert(all(all(edgesPos == sort(edgesPos,2))));
    assert(all(all(edgesNeg == sort(edgesNeg,2))));
    assert(all(all(edgesInBbox == sort(edgesInBbox,2))));
    assert(~any(ismember(edgesPos, edgesNeg, 'rows')));
    labels(ismember(edgesInBbox, edgesPos, 'rows')) = 1;
    labels(ismember(edgesInBbox, edgesNeg, 'rows')) = -1;
end
