function endings = buildEndings(axons)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    endings = table;
    
    endings.aggloId = repelem( ...
        reshape(1:numel(axons), [], 1), ...
        arrayfun(@(a) numel(a.endings), axons));
    endings.nodeId = cat(1, axons.endings);
end