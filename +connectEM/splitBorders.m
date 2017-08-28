function splitBorders(startidx, p)
    % load needed data
    temp = load(fullfile(p.saveFolder,'aggloState/axons_02.mat'));
    superagglos = temp.axons;
    graph = connectEM.loadAllSegmentationData(p);
    borderMeta = load([p.saveFolder 'globalBorder.mat'], 'borderCoM');

    % create bbox
    options.border = [2000; -2000];
    borderNm = repmat(options.border, 1, 3);
    borderVoxel = round(bsxfun(@times, 1./p.raw.voxelSize, borderNm));
    bboxSmall = p.bbox + borderVoxel';

    % create fake center of dataset field in borderCoM list for correspondences
    borderMeta.borderCoM(end+1,:)=round(mean(bboxSmall,2)');


    for idx = startidx : 50 : length(superagglos)
        nodesHere = superagglos(idx).nodes;
        segIdsHere = nodesHere(:, 4);
        flatten = @(x)x(:);
        edgesCandidates = sort([flatten(repelem(segIdsHere,cellfun('length',graph.neighbours(segIdsHere)))), cell2mat(graph.neighbours(segIdsHere))],2);
        edgesCandidatesProb = cell2mat(graph.neighProb(segIdsHere));
        edgesCandidatesBorderIdx = cell2mat(graph.neighBorderIdx(segIdsHere));
        only12 = @(x)x(:,1:2);
        edgeOverview = flipud(sortrows([double(edgesCandidates),edgesCandidatesProb,double(edgesCandidatesBorderIdx)],3));
        [~, filter] = unique(only12(edgeOverview),'rows');
        edgesCandidates=edgeOverview(filter,1:2);
        edgesCandidatesProb = edgeOverview(filter,3);
        edgesCandidatesBorderIdx = edgeOverview(filter,4);
        
        [~, usededges] = ismember(edgesCandidates,sort(reshape(segIdsHere(superagglos(idx).edges),[],2),2),'rows');
        % find whether edge is outside of bbox
        borderidxs = edgesCandidatesBorderIdx(usededges>0);
        borderidxs(isnan(borderidxs))=length(borderMeta.borderCoM); % correspondences
        borderPositions = double(borderMeta.borderCoM(borderidxs,:));
        outsideBbox = ~(all(bsxfun(@gt, borderPositions, bboxSmall(:, 1)'), 2) & ...
            all(bsxfun(@lt, borderPositions, bboxSmall(:, 2)'), 2));
        % remove all insulting edges
        edgesFinal = superagglos(idx).edges(usededges(usededges>0),:);
        if ~isempty(edgesFinal)
            edgesFinal(outsideBbox&edgesCandidatesProb(usededges>0)<0.98, :) = [];
        end
        C = Graph.findConnectedComponents(edgesFinal);
        superagglosBorderSplit{idx} = struct('nodes', {},'edges', {});
        for idx2 = 1 : length(C)
            superagglosBorderSplit{idx}(idx2).nodes = superagglos(idx).nodes(C{idx2},:);
            clear lookup
            lookup(C{idx2})=1:length(C{idx2});
            superagglosBorderSplit{idx}(idx2).edges = sort(reshape(lookup(edgesFinal(any(ismember(edgesFinal,C{idx2}),2),:)),[],2),2);
        end
        missingNodes = setdiff(1:size(superagglos(idx).nodes,1), cell2mat(C));
        superagglosBorderSplit{idx} = [superagglosBorderSplit{idx}, struct('edges',cell(1,length(missingNodes)),'nodes', arrayfun(@(x){superagglos(idx).nodes(x,:)},missingNodes))];
        calculateLength = @(x)max(pdist(bsxfun(@times, double(borderMeta.borderCoM(x, :)), p.raw.voxelSize)));
        filternan = @(x)x(~isnan(x));
        for idx2 = 1 : length(superagglosBorderSplit{idx})
            if length(superagglosBorderSplit{idx}(idx2).nodes(:,4)) < 5000
                axonLength{idx}(idx2) = max([-1, calculateLength(filternan(cell2mat(graph.neighBorderIdx(superagglosBorderSplit{idx}(idx2).nodes(:,4)))))]);
            else
                axonLength{idx}(idx2) = Inf;
            end
        end
        %comments = cellfun(@(x){cell(size(x,1), 1)},{superagglosBorderSplit{idx}.nodes});
        %writeNml(['a.nml'], writeSkeletonFromNodesAndEdges({superagglosBorderSplit{idx}.nodes}, {superagglosBorderSplit{idx}.edges},comments, repmat({'axon'},size(superagglosBorderSplit{idx})), repmat({[0 0 1 1]},size(superagglosBorderSplit{idx}))));
        %comments = cellfun(@(x){cell(size(x,1), 1)},{superagglos(idx).nodes});
        %writeNml(['b.nml'], writeSkeletonFromNodesAndEdges({superagglos(idx).nodes}, {superagglos(idx).edges},comments, repmat({'axon'},1,1), repmat({[0 0 1 1]},1,1)));

    end
    save([p.saveFolder 'splitBorder/splitBorder_' num2str(startidx)], 'superagglosBorderSplit','axonLength')
end
function dummy()
end