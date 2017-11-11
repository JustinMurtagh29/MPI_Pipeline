function out = prepareSplit(chiasmata, taskDefs, exits, taskIds, flights)
    % out = prepareSplit(chiasmata, taskDefs, exits, taskIds, flights)
    %   This function builds the data structure consumed by the chiasma
    %   splitting procedure in `connectEM.splitChiasmataMultiSuper`.
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    % sanity checks
    assert(size(taskDefs, 1) == size(exits, 1));
    assert(size(taskDefs, 1) == numel(taskIds));
    assert(all(ismember(flights.filenamesShort, taskIds)));
    
    out = struct;
    out.taskDef = data.taskDef;
    out.taskIds = taskIds;
    out.ff = flights;
    
    % build `queries` matrix
    out.queries = buildQueries(chiasmata, taskDefs, exits);
end

function queries = buildQueries(chiasmata, taskDefs, exits)
    taskCount = size(taskDefs, 1);
    queries = nan(taskCount, 1);
    
    for curIdx = 1:taskCount
        curExit = exits(curIdx, :);
        curPosition = taskDefs.position(curIdx, :) + 1;
        
        curChiasmata = chiasmata{curExit.aggloId};
        curCenterNodeId = curChiasmata.ccCenterIdx(curExit.chiasmaId);
        curExitNodeId = curChiasmata.queryIdx{curExit.chiasmaId};
        curExitNodeId = curExitNodeId(curExit.exitId);
        
        queries(curIdx, 1:3) = table2array(curExit);
        queries(curIdx, 4:6) = curPosition;
        queries(curIdx,   7) = curCenterNodeId;
        queries(curIdx,   8) = curExitNodeId;
    end
end