function cellToHdf5(file, group, cellArray)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    assert(iscell(cellArray));
    assert(isvector(cellArray));
    
    for curIdx = 1:numel(cellArray)
        curDset = fullfile(group, num2str(curIdx));
        arrayToHdf5(file, curDset, cellArray{curIdx});
    end
end
