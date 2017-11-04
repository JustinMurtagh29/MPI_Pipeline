function chiasmata = collectResults(chiasmaDir, chiasmaId, aggloCount)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    chiasmata = cell(aggloCount, 1);
    
    for curIdx = 1:aggloCount
        curFile = fullfile( ...
            chiasmaDir, ...
            sprintf('chiasmataX%s_%d', chiasmaId, floor(curIdx / 100)), ...
            sprintf('visX%s_%d', chiasmaId, curIdx), 'result.mat');

        curData = load(curFile);
        chiasmata{curIdx} = curData.output;
    end
end