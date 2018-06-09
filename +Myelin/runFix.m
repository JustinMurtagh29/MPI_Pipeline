function job = runFix(p, newPrefix)
    % job = runFix(p, newPrefix)
    %   Applies `Myelin.enforceMyelinSegments` to the entire dataset.
    %   For further information, check out the documentation of said
    %   function.
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    % TODO(amotta): This hack exploited the prefix of webKONSSOS cubes, and
    % is thus incompatible with the WKW format.
    validBackend = ...
       ~isfield(p.class, 'backend') || strcmp(p.class.backend, 'wkcube');
    assert(validBackend, [ ...
        'Cannot run %s with this storage backend. See also:\n', ...
        'https://gitlab.mpcdf.mpg.de/connectomics/pipeline/issues/37'], ...
        mfilename);
    
    taskInputArguments = arrayfun( ...
        @(local) {newPrefix, local.bboxSmall}, ...
        p.local(:), 'UniformOutput', false);
    
    job = Cluster.startJob( ...
        @Myelin.enforceMyelinSegments, ...
        taskInputArguments, ...
        'sharedInputs', {p}, ...
        'name', 'myelinFix', ...
        'cluster', { ...
            'memory', 12, ...
            'time', '12:00:00'});
end
