function prepareTrainingData(pT)
% Collects normalized training data to pT.initalGroundTruth and values needed for normalization of test data to pT.normValues

allIdx = [];
allLabels = [];
allWeights = [];
for i=1:length(pT.local)
	% Load data of current densly skeletonized training region
	[~, seg] = loadSegData(pT.local(i).segFile, pT.tileBorder);
	load(pT.local(i).edgeFile);
	load(pT.local(i).weightFile);
	skel.file = pT.local(i).trainFile;
	skel.bbox = pT.local(i).bboxSmall;
	% Extract labels for all segments that intersect with any of the dense skeletons
	[labelIdx, labels] = extractGroundTruthFromNml(seg, edges, skel);
	% Keep allWeights and just the labelled ones for normalization of test data & GP training respectively
	allWeights = [allWeights; weights];
	allIdx = [allIdx; labelIdx];
	allLabels = [allLabels; labels];
	clear seg edges weights skel labelIdx labels;
end

% Determine global 'whitening' values (0 and 1 of each feature correspond to min/max value in training set)
% Get min/max values (and precomputed value 'compFactor' for faster scaling to 0-1 on test data) of original data
minValues = min(allWeights,[],1);
maxValues = max(allWeights,[],1);
compFactor = 1./(maxValues-minValues);
save(pT.gp.normValues, 'minValues', 'maxValues', 'compFactor');

% Normalize allWeights 
allWeights = bsxfun(@minus,allWeights,minValues);
allWeights = bsxfun(@times,allWeights,compFactor);

% Split into training and test data (keep ~80% as training label)
nrTrainingSamples = round(length(allLabels)*.8);
nrTestSamples = length(allLabels) - nrTrainingSamples;
randIdx = [true(1,nrTrainingSamples) false(1,nrTestSamples)];
randIdx = randIdx(randperm(length(randIdx)));
allIdx = find(allIdx);
trainingData = allWeights(allIdx(randIdx),:);
trainingLabels = allLabels(randIdx);
testData = allWeights(allIdx(~randIdx),:);
testLabels = allLabels(~randIdx);
save(pT.gp.initalGroundTruth, 'trainingData', 'trainingLabels', 'testData', 'testLabels');

end

