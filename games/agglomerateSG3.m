function [agglo, probabilities, querried, merger, querriedEdges, mergerList] = agglomerateSG3(graph, com, seeds, segIds)
    % Agglomerates supervoxel given the graph, always highest probability for nrSteps

    % Initalization
    agglo = mat2cell(seeds, ones(length(seeds),1));
    probabilities = cell(size(agglo));
    querried = cell(size(agglo));
    merger = cell(size(agglo));
    querriedEdges = cell(size(agglo));
    mergerList = [];
    aggloFlag = true(size(agglo));

    while true 
        for s=1:length(seeds)
            if aggloFlag(s)
                % Remove merger from agglomerated segments
                thisAggloWithoutMerger = setdiff(agglo{s}, mergerList);
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
                if maxProb < .825
                    % Determine segment farest way from start com
                    startCoord = com(seeds(s),:);
                    aggloCoord = com(thisAggloWithoutMerger,:);
                    distances = sqrt(sum(bsxfun(@times, bsxfun(@minus, aggloCoord, startCoord), [11.24 11.24 28]).^2,2));
                    [~,maxDistIdx] = max(distances);
                    % Querry skeleton
                    addedId = querrySkeleton(graph, agglo{s}, segIds, thisAggloWithoutMerger(maxDistIdx));
                    querriedEdges{s}(end+1,1:2) = [thisAggloWithoutMerger(maxDistIdx) addedId];
                    querried{s}(end+1) = true;
                else
                    addedId = maxProbId;
                    querried{s}(end+1) = false;
                end
                % Add id to agglomeration
                agglo{s}(end+1) = addedId;
                probabilities{s}(end+1) = maxProb;
                % Next Heuristics: Exclude if 40 neighbours or more
                if length(graph.neighbours{addedId}) >= 40 
                    merger{s}(end+1) = true;
                    % Add last seed ID to ignore list (SegEM merger)
                    mergerList(end+1) = addedId;
                else
                    merger{s}(end+1) = false;
                end
            end
        end
        % New stop condition(s), stop growing of seeds as soon as two seeds meet
        for i=1:length(agglo)
            for j=i+1:length(agglo)
                overlap(i,j) = ~isempty(intersect(agglo{i},agglo{j}));
                if overlap(i,j)
                    aggloFlag(i) = false;
                    aggloFlag(j) = false;
                end
            end
        end
        % Stop agglomeration if all seeds have met
        if all(~aggloFlag)
            break;
        end       
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

