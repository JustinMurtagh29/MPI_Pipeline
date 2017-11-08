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
