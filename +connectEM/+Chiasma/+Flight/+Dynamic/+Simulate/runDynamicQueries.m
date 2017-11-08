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
chiasmaT.dynCorrect = false(size(chiasmaT, 1), 1);

chiasmaT.statLUT = cell(size(chiasmaT, 1), 1);
chiasmaT.statPart = cell(size(chiasmaT, 1), 1);

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
    curTracings = summaries(curAggloId).tracings{curChiasmaId};
    
    curStatLUT = curTracings.lut;
    curStatPart = reshape(curTracings.partition, 1, []);
    curStatPartValid = curTracings.partitionIsValid;
    
    %% build partition string for static solution
    curStatPartStr = arrayfun( ...
        @num2str, curStatPart, 'UniformOutput', false);
    curStatPartStr = strjoin(curStatPartStr, '-');
    
    curStatPartValidStr = {'invalid', 'valid'};
    curStatPartValidStr = curStatPartValidStr{1 + curStatPartValid};
    curStatPartStr = strcat(curStatPartStr, '-', curStatPartValidStr);
    
    %% collect outputs
    chiasmaT.dynLUT{curIdx} = curDynLUT;
    chiasmaT.statLUT{curIdx} = curStatLUT;
    chiasmaT.statPart{curIdx} = curStatPartStr;
end

chiasmaT.dynCorrect = cellfun( ...
    @isequal, chiasmaT.dynLUT, chiasmaT.statLUT);
chiasmaT.dynOverlaps = [];
clear cur*;

%% show results
chiasma4T = chiasmaT;
chiasma4T = chiasma4T(chiasma4T.nrExits == 4, :);

fprintf('\n');
fprintf('Total number of chiasmata: %d\n', size(chiasmaT, 1));
fprintf('Total number of chiasma exits: %d\n', sum(chiasmaT.nrExits));

fprintf('\n');
fprintf( ...
    'Total number of dynamic queries: %d (%.1f %%)\n', ...
    sum(chiasmaT.dynNrQueries), 100 * ...
    sum(chiasmaT.dynNrQueries) / sum(chiasmaT.nrExits));
fprintf( ...
    'Total number of dynamic query rounds: %d\n', ...
    max(chiasmaT.dynNrQueries));
fprintf( ...
    'Dynamic queries correct: %.1f %%\n', ...
    100 * mean(chiasmaT.dynCorrect));

temp = table;
[temp.nrExits, ~, tempUniRows] = unique(chiasmaT.nrExits);
temp.nrChiasmata = accumarray(tempUniRows, 1);
temp.nrCorrect = accumarray(tempUniRows, chiasmaT.dynCorrect);
temp.percCorrect = round(100 .* temp.nrCorrect ./ temp.nrChiasmata, 1);

fprintf('\n');
fprintf('Error rate vs. number of chiasma exits:\n\n');
disp(temp);

temp = table;
[temp.partition, ~ , tempUniRows] = unique(chiasma4T.statPart);
temp.nrChiasmata = accumarray(tempUniRows, 1);
temp.nrCorrect = accumarray(tempUniRows, chiasma4T.dynCorrect);
temp.percCorrect = round(100 .* temp.nrCorrect ./ temp.nrChiasmata, 1);

fprintf('\n');
fprintf('Error rate vs. chiasma partition (4-fold only):\n\n');
disp(temp);

temp = table;
[temp.partition, ~ , tempUniRows] = unique(chiasmaT.statPart);
temp.nrChiasmata = accumarray(tempUniRows, 1);
temp.nrCorrect = accumarray(tempUniRows, chiasmaT.dynCorrect);
temp.percCorrect = round(100 .* temp.nrCorrect ./ temp.nrChiasmata, 1);

fprintf('\n');
fprintf('Error rate vs. chiasma partition:\n\n');
disp(temp);

clear chiasma4T;
clear temp*;