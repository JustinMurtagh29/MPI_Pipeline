function nullMaps = nullDensityMaps(asiAreas, varargin)
    opts = struct;
    opts.numMaps = 1;
    opts.scheduler = '';
    opts.asiGroups = [];
    opts = Util.modifyStruct(opts, varargin{:});
    
    % NOTE(amotta): Make sure that user has specified bandwidth for kernel
    % density estimate. Otherwise we're comparing apples and oranges...
    assert(isfield(opts, 'bandWidth'));
    assert(not(isempty(opts.bandWidth)));
    
    asiGroups = {1:numel(asiAreas)};
    if ~isempty(opts.asiGroups)
        assert(isequal(numel(asiAreas), numel(opts.asiGroups)));
       [~, ~, asiGroups] = unique(opts.asiGroups);
       
        asiGroups = accumarray( ...
            reshape(asiGroups, [], 1), ...
            reshape(1:numel(asiAreas), [], 1), ...
            [], @(ids) {ids(:)});
    end
    
    args = arrayfun(@(i) {i}, 1:opts.numMaps, 'UniformOutput', false);
    
    if isempty(opts.scheduler)
        nullMaps = cellfun( ...
            @(args) nullDensityMap( ...
                asiAreas, asiGroups, args{:}, varargin{:}), ...
            args, 'UniformOutput', false);
    else
        job = Cluster.startJob( ...
            @nullDensityMap, args, ...
            'numOutputs', 1, ...
            'name', mfilename, ...
            'taskGroupSize', 20, ...
            'sharedInputs', [{asiAreas, asiGroups}, varargin], ...
            'sharedInputsLocation', [1:2, 3 + (1:numel(varargin))], ...
            'cluster', { ...
                'memory', 12, ...
                'priority', 100, ...
                'scheduler', opts.scheduler});
        Cluster.waitForJob(job);

        nullMaps = fetchOutputs(job);
        nullMaps = cat(1, nullMaps{:});
    end
    
    nullMaps = reshape(nullMaps, 1, 1, []);
    nullMaps = cell2mat(nullMaps);
end

function nullMap = nullDensityMap(asiAreas, asiGroups, rngSeed, varargin)
    rng(rngSeed);
    
    asiAreaPairs = zeros(0, 2);
    for curIdx = 1:numel(asiGroups)
        curIds = asiGroups{curIdx};
        
        curPairs = 2 * floor(numel(curIds) / 2);
        curPairs = randperm(numel(curIds), curPairs);
        curPairs = curIds(reshape(curPairs, [], 2));
        
        asiAreaPairs = cat(1, asiAreaPairs, curPairs);
    end
    
    asiAreaPairs = asiAreas(asiAreaPairs);
    nullMap = connectEM.Consistency.densityMap(asiAreaPairs, varargin{:});
end
