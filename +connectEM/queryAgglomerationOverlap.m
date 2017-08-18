function [partition, queryOverlap] = queryAgglomerationOverlap(agglos, segmentsLeftover, uniqueSegments, neighboursStartNode)
    if isfield(agglos,'nodes') % new representation
        agglos = cellfun(@(x) x(:,4),{agglos.nodes},'uni',0);
        if length(agglos)>size(agglos,1)
            agglos = agglos';
        end
    end
    if isfield(segmentsLeftover,'nodes') % new representation
        segmentsLeftover = cellfun(@(x) x(:,4),{segmentsLeftover.nodes},'uni',0);
    end
    
    % Find start agglomeration and overlap of each query 

    % First generate different data structure for speedup of loop below (4.6s per iteration with method as in debugQueryA...)
    % Also add all leftover segments as single equivlaence classes to also collect them if sufficent information
    partition.segIds = cat(1, agglos{:}, segmentsLeftover);
    partition.eqClass = (1:length(partition.segIds))';
    sizeOfComponent = [0; cumsum(cellfun(@length, agglos))];
    for i=1:length(sizeOfComponent)-1
        partition.eqClass(sizeOfComponent(i)+1:sizeOfComponent(i+1)) = i;
    end
    clear sizeOfComponent;

    queryOverlap.ends = cellfun(@(x)determineGroupedAgglomerationOverlap(partition, x), uniqueSegments);
    queryOverlap.start = cellfun(@(x)determineGroupedAgglomerationOverlap(partition, x), neighboursStartNode);

end

function result = determineGroupedAgglomerationOverlap(partition, segmentStruct)
    persistent index
    if isempty(index)
        index = 0;
    else
        index =  index+1;
    end
    disp(index);
    % Determine agglomeration that each segment belongs to and evidence (occurence) for each
    [Lia, Locb] = ismember(partition.segIds, segmentStruct.segIds);
    eqClass = partition.eqClass(Lia)';
    occurences = segmentStruct.occurences(Locb(Lia))';
    % eqClass2 = [];
    % occurences2 = [];
    % for i=1:length(segmentStruct.segIds)
    %     idx = partition.segIds == segmentStruct.segIds(i);
    %     if any(idx)
    %         eqClass2(end+1) = partition.eqClass(idx);
    %         occurences2(end+1) = segmentStruct.occurences(i);
    %     end
    % end
    % assert(isequal(sortrows([eqClass' occurences']),sortrows([eqClass2' occurences2'])));
    result.eqClasses = unique(eqClass);
    if ~isempty(result.eqClasses)
        for i=1:length(result.eqClasses)
            result.occurences(i) = sum(occurences(eqClass == result.eqClasses(i)));
        end
        result = sortAccordingToOccurence(result);
    else
        result.occurences = [];
    end
end

function in = sortAccordingToOccurence(in)
    [in.occurences, idxPerm] = sort(in.occurences, 'descend');
    in.eqClasses = in.eqClasses(idxPerm);
end

