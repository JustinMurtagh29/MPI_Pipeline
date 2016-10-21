function job = runMyelinFix(p, newPrefix)
    % job = runMyelinFix(p, newPrefix)
    %   Applies `runMyelinFixBox` to the entire dataset. For
    %   further information, check out the documentation of
    %   said function.
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    xIds = p.bbox(1, 1):128:p.bbox(1, 2);
    yIds = p.bbox(2, 1):128:p.bbox(2, 2);
    zIds = p.bbox(3, 1):128:p.bbox(3, 2);
    
    % do *NOT* change X and Y!
    [xMinVec, yMinVec, zMinVec] = meshgrid(yIds, xIds, zIds);
    
    % build job inputs
    boxCount = numel(xMinVec);
    buildBox = @(idx) [ ...
        xMinVec(idx), xMinVec(idx) + 127;
        yMinVec(idx), yMinVec(idx) + 127;
        zMinVec(idx), zMinVec(idx) + 127];
    
    taskInputArguments = arrayfun( ...
        @(idx) {{p, newPrefix, buildBox(idx)}}, 1:boxCount);
    job = startCPU(@runMyelinFixBox, taskInputArguments, 'myelinFix');
end
