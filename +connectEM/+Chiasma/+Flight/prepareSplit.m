function out = prepareSplit(chiasmata, taskDef, exits, taskIds, flights)
    % out = prepareSplit(chiasmata, taskDef, exits, taskIds, flights)
    %   This function builds the data structure consumed by the chiasma
    %   splitting procedure in `connectEM.splitChiasmataMultiSuper`.
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    % sanity checks
    assert(size(taskDef, 1) == size(exits, 1));
    assert(size(taskDef, 1) == numel(taskIds));
    assert(all(ismember(flights.filenamesShort, taskIds)));
    
    out = struct;
    out.ff = flights;
    out.taskDef = taskDef;
    out.taskIds = taskIds;
    
    % build `queries` matrix
    out.queries = buildQueries(chiasmata, taskDef, exits);
end

function queries = buildQueries(chiasmata, taskDef, exits)
    taskCount = size(taskDef, 1);
    queries = nan(taskCount, 1);
    
    for curIdx = 1:taskCount
        curExit = exits(curIdx, :);
        curPosition = taskDef.position(curIdx, :) + 1;
        
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