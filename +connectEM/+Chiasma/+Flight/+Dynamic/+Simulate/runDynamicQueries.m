% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>

%% run chiasma splitting on fully queried chiasmata for ground truth
summaries = connectEM.Chiasma.Flight.Dynamic.Simulate.runStaticQueries();

%% build fake chiasma table
knownNrExits = reshape({summaries.nrExits}, [], 1);
knownNrChiasmata = cellfun(@numel, knownNrExits);
knownAggloCount = numel(knownNrChiasmata);

chiasmaT = table;
chiasmaT.aggloId = repelem( ...
    reshape(1:knownAggloCount, [], 1), knownNrChiasmata);
chiasmaT.chiasmaId = cell2mat(arrayfun( ...
    @(count) reshape(1:count, [], 1), ...
    knownNrChiasmata, 'UniformOutput', false));
chiasmaT.nrExits = cell2mat(knownNrExits);

%% prepare overlap structure
overlaps = ...
    cellfun( ...
        @(allNrExits) arrayfun( ...
            @(nrExits) zeros(nrExits, 2), ...
            allNrExits, 'UniformOutput', false), ...
        knownNrExits, 'UniformOutput', false);

clear known*;

%% simulate querying
curRoundIdx = 1;

while true
    curExits = connectEM.Chiasma.Flight.selectExits(chiasmaT, overlaps, 1);
    if isempty(curExits); break; end
    
    fprintf('Query round %d\n', curRoundIdx);
    fprintf('  Number of queries: %d\n', size(curExits, 1));
    
    fprintf('  Simulating queries... ');
    for curIdx = 1:size(curExits, 1)
        curAggloId = curExits.aggloId(curIdx);
        curChiasmaId = curExits.chiasmaId(curIdx);
        curExitId = curExits.exitId(curIdx);

        % simulate query results
        curOverlap = summaries(curAggloId).tracings{curChiasmaId}.overlaps{curExitId};
        overlaps{curAggloId}{curChiasmaId}(curExitId, :) = curOverlap;
    end
    fprintf('done!\n');
    
    curRoundIdx = curRoundIdx + 1;
end

clear cur*;

%% collect overlaps and build LUT
chiasmaT.dynOverlaps = cat(1, overlaps{:});
chiasmaT.dynNrQueries = cellfun( ...
    @(o) sum(o(:, 1) > 0), chiasmaT.dynOverlaps);

chiasmaT.dynLUT = cell(size(chiasmaT, 1), 1);
chiasmaT.statLUT = cell(size(chiasmaT, 1), 1);

for curIdx = 1:size(chiasmaT, 1)
    curOverlaps = chiasmaT.dynOverlaps{curIdx};
    curNrExits = size(curOverlaps, 1);
    
    %% build lut
    % This code was copied from `connectEM.splitChiasmataMultiLogic` at
    % git commit bf36ada78e62970be738fd6ee2351e7573f63740.
    
    % build effective edges
    curEdges = sort(curOverlaps, 2);
    curEdges(~all(curEdges, 2), :) = [];
    curEdges = unique(curEdges, 'rows');
    curEdges = reshape(curEdges, [], 2);
    
    % find parition
    curAdjMat = sparse( ...
        curEdges(:, 2), curEdges(:, 1), ...
        true, curNrExits, curNrExits);
   [~, curDynLUT] = graphconncomp(curAdjMat, 'Directed', false);
    curDynLUT = reshape(curDynLUT, [], 1);
    
    %% retrieve LUT from static querying
    curAggloId = chiasmaT.aggloId(curIdx);
    curChiasmaId = chiasmaT.chiasmaId(curIdx);
    curStatLUT = summaries(curAggloId).tracings{curChiasmaId}.lut;
    
    %% collect outputs
    chiasmaT.dynLUT{curIdx} = curDynLUT;
    chiasmaT.statLUT{curIdx} = curStatLUT;
end

clear cur*;