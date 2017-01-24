function job = runFix(p, newPrefix)
    % job = runFix(p, newPrefix)
    %   Applies `+Myelin/runFixInBox` to the entire dataset.
    %   For further information, check out the documentation
    %   of said function.
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    boxCount = numel(p.local);
    taskInputArguments = arrayfun( ...
        @(idx) {{p, newPrefix, p.local(idx).bboxSmall}}, 1:boxCount);
    job = startCPU(@Myelin.enforceMyelinSegments, taskInputArguments, 'myelinFix');
end

