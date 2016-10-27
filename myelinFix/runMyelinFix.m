function job = runMyelinFix(p, newPrefix)
    % job = runMyelinFix(p, newPrefix)
    %   Applies `runMyelinFixBox` to the entire dataset. For
    %   further information, check out the documentation of
    %   said function.
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    boxCount = numel(p.local);
    taskInputArguments = arrayfun( ...
        @(idx) {{p, newPrefix, p.local(idx).bboxSmall}}, 1:boxCount);
    
    job = startCPU(@runMyelinFixBox, taskInputArguments, 'myelinFix');
end
