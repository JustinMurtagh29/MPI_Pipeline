param = load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat');
param = param.p;
dataDir = fullfile(param.saveFolder, 'aggloState');

superAgglos = load(fullfile(dataDir, 'axons_04.mat'));
origAgglos = arrayfun(@Agglo.fromSuperAgglo, superAgglos.axons, 'uni', 0);
superAgglos = superAgglos.axons(superAgglos.indBigAxons);

endings = load(fullfile(dataDir, 'axonEndings.mat'));
endingOverlap = load(fullfile(dataDir, 'axonEndingOverlaps.mat'));

clusterSizes = cellfun(@max, endings.borderClusters);

% Which endings are linked together
linkages = cellfun(@(x,y)edgeCreator(x,y), endingOverlap.startEndingOverlaps, endingOverlap.endEndingOverlaps, 'uni', 0);
linkagesFlat = cat(1, linkages{:});

% Which endings belong to which agglo
clusterSizes = cellfun(@max, endings.borderClusters);
clusterLookup = repelem(1:numel(clusterSizes), clusterSizes);

caseDistinctions = zeros(size(linkagesFlat,1),1);
% Choose 1 for no attachment
caseDistinctions(isnan(linkagesFlat(:,1)) & isnan(linkagesFlat(:,2))) = 1;
% Choose 2 for attachment at end but not at start
caseDistinctions(isnan(linkagesFlat(:,1)) & ~isnan(linkagesFlat(:,2))) = 2;
% Choose 3 for attachment to start but not at end
caseDistinctions(~isnan(linkagesFlat(:,1)) & isnan(linkagesFlat(:,2))) = 3;
% Find indices did match criteria so far
idxMatched = any(isnan(linkagesFlat),2);
% Choose 4 for attachted to start not at ending and to end not at ending
caseDistinctions(linkagesFlat(:,1) == 0 & linkagesFlat(:,2) == 0 & ~idxMatched) = 4;
% Choose 5 for attachted to start at ending and to end not at ending
caseDistinctions(linkagesFlat(:,1) ~= 0 & linkagesFlat(:,2) == 0 & ~idxMatched) = 5;
% Choose 6 for attachted to start not at ending and to end at ending
caseDistinctions(linkagesFlat(:,1) == 0 & linkagesFlat(:,2) ~= 0 & ~idxMatched) = 6;
% Find indices did not match criteria so far
idxMatched = idxMatched | any(linkagesFlat == 0,2);
% Choose 7 for self attachment (within ending)
caseDistinctions((linkagesFlat(:,1) == linkagesFlat(:,2)) & ~idxMatched) = 7;
% Find indices did not match criteria so far
idxMatched = idxMatched | (linkagesFlat(:,1) == linkagesFlat(:,2));
% Choose 8 for self attachment (within agglo)
caseDistinctions((clusterLookup(linkagesFlat(:,1)) == clusterLookup(linkagesFlat(:,2))) & ~idxMatched) = 8;
% Find indices did not match criteria so far
idxMatched = idxMatched | (clusterLookup(linkagesFlat(:,1)) == clusterLookup(linkagesFlat(:,2)));
% Choose 9 for normal attachment
caseDistinctions(~idxMatched) = 9;

% Output count & frequency of each case
tabulate(caseDistinctions)

