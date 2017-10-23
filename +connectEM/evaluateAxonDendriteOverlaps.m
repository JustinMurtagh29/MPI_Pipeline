% This function takes a bunch of axon super-agglomerates and checks if the
% flight paths contained therein have significant overlaps with some of the
% dendrites agglomerates.
%
% The hope is that this helps us to identify axon fragments which ended up
% in the wrong bucket. It could be that these fragments are a major source
% of "dendrite" endings. Identifying these fragments could potentially help
% us to get rid of unnecessary queries.
%
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>

clear;
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
axonFile = fullfile(rootDir, 'chiasmataSplitting/20171009T193744-kmb-on-axons-6c/outputs/20171023T115300_results.mat');
dendFile = fullfile(rootDir, 'aggloState/dendrites_03_v2_splitmerged.mat');

%% load input data
param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

axons = load(axonFile);
axonIds = find(axons.indBigAxons);
axons = axons.axons(axonIds);

%% reconstruct flight paths
axonCount = numel(axonIds);
axonFlights = cell(axonCount, 1);
axonFlights(:) = {zeros(0, 6)};

for curIdx = 1:axonCount
    curAxon = axons(curIdx);
    
    % work around empty edges
    curAxon.edges = reshape(curAxon.edges, [], 2);
    
    % extract flight nodes and edges
    curFlightNodes = find(isnan(curAxon.nodes(:, 4)));
   [~, curFlightEdges] = ismember(curAxon.edges, curFlightNodes);
    curFlightEdges = curFlightEdges(all(curFlightEdges, 2), :);
    curFlightEdges = sort(curFlightEdges, 2);
    
    % group into paths
    curAdj = sparse( ...
        curFlightEdges(:, 2), curFlightEdges(:, 1), ...
        true, numel(curFlightNodes), numel(curFlightNodes));
   [curFlightCount, curLUT] = graphconncomp(curAdj, 'Directed', false);
   
    axonFlights{curIdx} = horzcat( ...
       repelem(curIdx, numel(curFlightNodes), 1), ...
       curFlightNodes(:), curLUT(:), curAxon.nodes(curFlightNodes, 1:3));
end

axonFlights = cell2mat(axonFlights);
axonFlights = transpose(axonFlights);

%% add 26-neighbours
posOff = cell(3, 1);
[posOff{:}] = ndgrid(-1:1, -1:1, -1:1);
posOff = cellfun(@(o) reshape(o, 1, []), posOff, 'UniformOutput', false);
posOff = cell2mat(posOff(:));

posOff = cat(1, zeros(size(axonFlights, 1) - 3, size(posOff, 2)), posOff);
axonFlights = reshape(axonFlights, 6, 1, []);
axonFlights = bsxfun(@plus, axonFlights, posOff);

axonFlights = reshape(axonFlights, 6, []);
axonFlights = transpose(axonFlights);

%% look up segment IDs
axonFlightSegIds = Seg.Global.getSegIds( ...
    param, axonFlights(:, (end - 2):end));