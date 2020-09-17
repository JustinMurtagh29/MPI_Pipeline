function [flights] = loadFromRoboTDefs(tdefOutputDir, tdefFNames, roboRunIds, param)
%% Written by
%    Martin Schmidt <martin.schmidt@brain.mpg.de>

% Read taskdefs and tracings:
tdefs = roboUtil.readBatchedTDefs(tdefOutputDir, tdefFNames, roboRunIds);

% Get main tracing (i.e. without validation path) of validated tracings only:
tracings = cellfun( ...
    @(tdef) tdef.getValidatedTracings(0), ...
    tdefs, 'uni', false);
tripods = cellfun( ...
    @(tracing) tracing.tripod, ...
    tracings, 'uni', false);
tripods = cat(1, tripods{:});
mainTracingTripods = cellfun( ...
    @(tripod) tripod{1}, ...
    tripods, 'uni', false);

% Get coordinates (i.e. 'nodes')
nodes = cell(numel(mainTracingTripods), 1);
nonEmptyTripods = cellfun(@numel, mainTracingTripods) > 0;
nodes(nonEmptyTripods) = cellfun( ...
    @(tripod) ...
        round(squeeze(tripod(:, 1, :)) ./ ...
        tdefs{1}.dset_params.scale) + [1, 1, 1], ... Matlab offset [1, 1, 1]
    mainTracingTripods(nonEmptyTripods), 'uni', false);


% Downsample to get reasonably spaced (not too dense) coordinates
nodes(nonEmptyTripods) = cellfun(@(curNodes) taskdef.dsample(curNodes, 10), nodes(nonEmptyTripods), 'uni', false);
nodes = cellfun(@double, nodes, 'uni', false);

% Get correct sorting using the global Task Id:
globTaskIds = Util.cellVertCat(cellfun(@(tdef) tdef.tasks.identifier.globTaskId, tdefs, 'uni', false));
[globTaskIdsSrtd, idcsSrtd] = sort(globTaskIds);

% Construct flights struct
flights = struct();
flights.nodes = nodes(idcsSrtd);
flights.startNode = flights.nodes;
nonEmptyFlights = cellfun(@numel, flights.startNode) > 0;
flights.startNode(nonEmptyFlights) = cellfun(@(x) x(1, :), flights.startNode(nonEmptyFlights), 'uni', false);

[flights.segIds, flights.neighbours] = connectEM.lookupSkelGT(param, flights);

% Add segIds, neighbours, filenames, filenamesShort, startNode, comments, globTaskId
numFlights = numel(nodes);
flights.comments = cell(numFlights, 1);
flights.filenames = cell(numFlights, 1);
flights.filenamesShort = globTaskIdsSrtd; % > dummy to be reassociated with globTaskId
flights.globTaskIds = globTaskIdsSrtd;

end
