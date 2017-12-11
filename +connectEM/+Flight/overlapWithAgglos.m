function overlaps = overlapWithAgglos(param, flights, agglos, varargin)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    %% configuration
    config = struct;
    config.minStartEvidence = 13;
    config.minEndEvidence = 2 * 27;
    
    config = Util.modifyStruct(config, varargin{:});
    
    %% determine overlap
    maxSegId = Seg.Global.getMaxSegId(param);
    aggloLUT = Agglo.buildLUT(maxSegId, agglos);
    aggloLUT = vertcat(0, reshape(aggloLUT, [], 1));
    
    flightCount = numel(flights.filenames);
    overlaps = cell(flightCount, 2);

    for curIdx = 1:flightCount
        curSegIds = 1 + horzcat( ...
            flights.segIds{curIdx}, ...
            flights.neighbours{curIdx});
        curEndId = aggloLUT(curSegIds);

        % find start agglo (> 1/2 node evidence)
       [curStartAgglos, ~, curStartEvidence] = unique(curEndId(1, :));
        curStartEvidence = accumarray(curStartEvidence, 1);
        curStartAgglos = curStartAgglos( ...
            curStartEvidence >= config.minStartEvidence);
        curStartAgglos = setdiff(curStartAgglos, 0);

        % find end agglos
       [curEndAgglos, ~, curEndEvidence] = unique(curEndId(:));
        curEndEvidence = accumarray(curEndEvidence, 1);
        curEndAgglos = curEndAgglos( ...
            curEndEvidence >= config.minEndEvidence);
        curEndAgglos = setdiff(curEndAgglos, union(0, curStartAgglos));

        overlaps{curIdx, 1} = curStartAgglos;
        overlaps{curIdx, 2} = curEndAgglos;
    end
end