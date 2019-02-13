function nullMaps = nullDensityMaps(asiAreas, varargin)
    opts = struct;
    opts.numMaps = 1;
    opts = Util.modifyStruct(opts, varargin{:});
    
    % NOTE(amotta): Make sure that user has specified bandwidth for kernel
    % density estimate. Otherwise we're comparing apples and oranges...
    assert(isfield(opts, 'bandWidth'));
    assert(not(isempty(opts.bandWidth)));
    
    args = arrayfun(@(i) {i}, 1:opts.numMaps, 'UniformOutput', false);
    sharedArgs = {asiAreas, varargin{:}}; %#ok
    
    job = Cluster.startJob( ...
        @nullDensityMap, args, ...
        'name', mfilename, ...
        'sharedInputs', sharedArgs, ...
        'sharedInputsLocation', [1, 2 + (1:numel(varargin))], ...
        'cluster', {'priority', 100, 'memory', 12}, ...
        'numOutputs', 1, 'taskGroupSize', 20);
    Cluster.waitForJob(job);
    
    nullMaps = fetchOutputs(job);
    nullMaps = cat(1, nullMaps{:});
    
    nullMaps = reshape(nullMaps, 1, 1, []);
    nullMaps = cell2mat(nullMaps);
end

function nullMap = nullDensityMap(asiAreas, rngSeed, varargin)
    rng(rngSeed);
    
    asiAreaPairs = 2 * floor(numel(asiAreas) / 2);
    asiAreaPairs = randperm(numel(asiAreas), asiAreaPairs);
    asiAreaPairs = reshape(asiAreas(asiAreaPairs), [], 2);
    
    nullMap = connectEM.Consistency.densityMap(asiAreaPairs, varargin{:});
end
