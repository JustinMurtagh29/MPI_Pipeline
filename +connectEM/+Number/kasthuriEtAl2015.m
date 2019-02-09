% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
info = Util.runInfo();
Util.showRunInfo(info);

% See page 653, Cell 162, 648-661, July 2015:
% "The spreadsheet shows that the connectivity is highly skewed towards
% excitatory elements: [..], and 93% (1,308/1,407) of the axons are
% excitatory."
axonCount = 1407;

%% Load synapses
xlsFile = fileparts(mfilename('fullpath'));
xlsFile = fullfile(xlsFile, 'kasthuriEtAl2015CellTableS1.xls');

data = xlsread(xlsFile);
data = data(3:end, :); % Get rid of header

colToId = @(c) 1 + double(lower(c)) - double('a');

synT = table;
synT.id = data(:, colToId('A'));
synT.axonId = data(:, colToId('I'));
synT.axonType = categorical( ...
    data(:, colToId('K')), [0, 1, 2, -1], ...
    {'Excitatory', 'Inhibitory', 'Myelinated', 'Unknown'});
synT.dendId = data(:, colToId('J'));
synT.ontoSpine = logical(data(:, colToId('R')));

% Sanity check
assert(isequal(synT.id, reshape(1:height(synT), [], 1)));

%% Calculate average number of synapses
[uniAxonIds, ~, synT.relAxonId] = unique(synT.axonId);

% Sanity check
assert(numel(uniAxonIds) <= axonCount);

axonT = table;
axonT.id = nan(axonCount, 1);
axonT.id(1:numel(uniAxonIds)) = uniAxonIds;
axonT.relId = reshape(1:axonCount, [], 1);
axonT.synCount = accumarray(synT.relAxonId, 1, size(axonT.id));

%% Evaluation
numSynsPerAxonMean = mean(axonT.synCount) %#ok
numSynsPerAxonStd = std(axonT.synCount) %#ok
numAxonsWithAtLeastTenSynapses = sum(axonT.synCount >= 10) %#ok

% Count multi-synaptic connections
[~, ~, numMultiSynPairs] = unique(synT(:, {'axonId', 'dendId'}), 'rows');
numMultiSynConnections = sum(accumarray(numMultiSynPairs, 1) > 1) %#ok

% Count multi-synaptic exc. connections onto spine heads
numExcShMultiSynPairs = synT( ...
    synT.axonType == 'Excitatory' ...
  & synT.ontoSpine, {'axonId', 'dendId'});
[~, ~, numExcShMultiSynPairs] = unique(numExcShMultiSynPairs, 'rows');
numExcShMultiSynConnections = sum(accumarray(numExcShMultiSynPairs, 1) > 1) %#ok
