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

% I thought it was weird that we have no self attachment (case 7 or 8), so I investigated:
% startAgglo already excluded from endAgglo in an earlier function
m = load(fullfile(dataDir, 'axonPostQueryAnalysisState.mat'));
startAgglo = arrayfun(@(x)x.eqClasses(x.occurences > 13), m.queryOverlap.start, 'uni', 0);
endAgglo = arrayfun(@(x)x.eqClasses(x.occurences > 53), m.queryOverlap.ends, 'uni', 0);
sum(cellfun(@(x,y)numel(intersect(x,y)), endAgglo, startAgglo, 'uni', 0));

% Check that border positions of new tasks make sense (needs linearized version of thisborderPositions from generateAxonQueries)

% Determine all endings hit by at least one flight path
endingOverlap = load(fullfile(dataDir, 'axonEndingOverlaps.mat'));
startEndings = unique(cell2mat(endingOverlap.startEndingOverlaps));
endEndings = unique(cell2mat(endingOverlap.endEndingOverlaps));
totalEndings = union(startEndings,endEndings);

% Determine open endings
endings = load(fullfile(dataDir, 'axonEndings.mat'));
openEndings = setdiff(1:sum(cellfun(@max, endings.borderClusters)), totalEndings);
numel(openEndings)

% Create lookup which agglo and which cluster within agglo an ending belongs to
clusterSizes = cellfun(@max, endings.borderClusters);
clusterLookup = repelem(1:numel(clusterSizes), clusterSizes);
withinClusterLookup = arrayfun(@(x)(1:x)', clusterSizes, 'uni', 0);
withinClusterLookup = cat(1, withinClusterLookup{:});

% Get positions of queried locations in cluster
theseBorderPositions = cellfun(@(x)cat(1, x{:}), thisborderPositions, 'uni', 0);
theseBorderPositions = cat(1, theseBorderPositions{:});

% Get positions of all members of each ending cluster not hit
positions = arrayfun(@(x)endings.borderPositions{clusterLookup(x)}...
    (endings.borderClusters{clusterLookup(x)} == withinClusterLookup(x),:), ...
    openEndings, 'uni', 0);
clusterSizes = cellfun(@(x)size(x,1), positions);
positions = cat(1, positions{:});

% Define exclusion zone near dataset center
options.border = [3000; -3000];
borderNm = repmat(options.border, 1, 3);
borderVoxel = round(bsxfun(@times, 1./param.raw.voxelSize, borderNm));
bboxSmall = param.bbox + borderVoxel';

% Exclude positions outside
idxWithin = all(bsxfun(@gt, positions, bboxSmall(:, 1)'),2) ...
    & all(bsxfun(@lt, positions, bboxSmall(:, 2)'),2);
positions = positions(idxWithin,:);

% Test whether positions are the same
suprise = setdiff(theseBorderPositions, positions, 'rows');
% None found, hence all seems fine

