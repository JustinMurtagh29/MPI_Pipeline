function [agglo, probabilities, stats, querriedEdges, mergerList] = agglomerateSG3(graph, com, seeds, termNodes, segIds)
    % Agglomerates supervoxel given the graph, always highest probability for nrSteps

    % Initalization
    display('Initalization');
    tic;
    agglo = seeds;
    probabilities = cell(size(seeds));
    querriedEdges = cell(size(seeds));
    mergerList = cell(size(seeds));
    stats = struct('querried', {}, 'merger', {});
    stats(4).querried = []; 
    toc;

    display('Agglomerating supervoxel');
    tic;
    for s=1:length(seeds)
        while isempty(intersect(agglo{s}, termNodes{s})) & length(agglo{s}) < 400
            % Remove merger from agglomerated segments
            thisAggloWithoutMerger = setdiff(agglo{s}, mergerList{s});
            % Determine neighbours of all agglomerated objects
            theseNeighbours = [graph.neighbours{thisAggloWithoutMerger}];
            % Determine probability of all neighbours
            theseProbabilities = cat(1,graph.neighProb{thisAggloWithoutMerger});
            % Remove self edges (within component)
            selfEdgeIdx = ismember(theseNeighbours, agglo{s});                      
            theseNeighbours(selfEdgeIdx) = [];
            theseProbabilities(selfEdgeIdx) = [];
            % Determine maxium probability and ID to agglomerated segments
            [maxProb, maxProbIdx] = max(theseProbabilities);
            % ID about to be agglomerated
            maxProbId = theseNeighbours(maxProbIdx);
            % Heuristics: Querry skeleton if sth sth sth dark side
            if maxProb < .85
                % Querry skeleton for continuation at point farest way from start com
                startCoord = com(seeds{s}(1),:);
                aggloCoord = com(agglo{s},:);
                distances = sqrt(sum(bsxfun(@times, bsxfun(@minus, aggloCoord, startCoord), [11.24 11.24 28]).^2,2));
                [~,maxDistIdx] = max(distances);
                addedId = querrySkeleton(graph, agglo{s}, segIds{s}, agglo{s}(maxDistIdx));
                querriedEdges{s}(end+1,1:2) = [agglo{s}(maxDistIdx) addedId];
                stats(s).querried(end+1) = true;
            else
                addedId = maxProbId;
                stats(s).querried(end+1) = false;
            end
            % Add id to agglomeration
            agglo{s}(end+1) = addedId;
            probabilities{s}(end+1) = maxProb;
            % Next Heuristics: Exclude if 40 neighbours or more
            if length(graph.neighbours{addedId}) >= 30 
                stats(s).merger(end+1) = true;
                % Add last seed ID to ignore list (SegEM merger)
                mergerList{s}(end+1) = addedId;
            else
                stats(s).merger(end+1) = false;
            end
        end
        % Display progress
        Util.progressBar(s, length(seeds));
    end

end

function newId = querrySkeleton(graph, agglo, segIds, querryId)
    % Find segment ID hit by skeleton closest to agglomeration
    tempIds = setdiff(segIds, agglo);
    found = false;
    currentIds = querryId;
    currentProb = 1;
    while ~found
        ids = [graph.neighbours{currentIds}];
        prob = [];
        for i=1:length(currentIds)
            prob{i} = graph.neighProb{currentIds(i)}'.*currentProb(i);
        end
        prob = cell2mat(prob);
        idx = ismember(ids, tempIds);
        if any(idx)
            found = true;
            idx = find(idx);
            [~,maxIdx] = max(prob(idx));
            newId = ids(idx(maxIdx));
        else
            [currentIds, idx2] = setdiff(ids, agglo);
            currentProb = prob(idx2);
        end
    end
end

