function prepareTrainingData(pT)
% Collects normalized training data to pT.initalGroundTruth and values needed for normalization of test data to pT.normValues

allIdx = [];
allLabels = [];
allWeights = [];
for i=1:length(pT.local)
	% Load data of current densly skeletonized training region
	load(pT.local(i).segFile);
	seg = seg(1-pT.tileBorder(1,1):end-pT.tileBorder(1,2),...
		1-pT.tileBorder(2,1):end-pT.tileBorder(2,2),...
		1-pT.tileBorder(3,1):end-pT.tileBorder(3,2));
	load(pT.edgeFile);
	load(pT.weightFile);
	skel.file = pT.local(i).trainFile;
	skel.bbox = pT.local(i).bboxSmall;
	% Extract labels for all segments that intersect with any of the dense skeletons
	[labelIdx, labels] = extractGroundTruthFromNml(seg, edges, skel);
	% Keep allWeights and just the labelled ones for normalization of test data & GP training respectively
	allWeights = [allWeights weights];
	allIdx = [allIdx labelIdx];
	allLabels = [allLabels labels];
	clear seg edges weights skel labelIdx labels;
end

%% Determine global 'whitening' values (0 and 1 of each feature correspond to min/max value in training set)
% Get min/max values (and precomputed value 'compFactor' for faster scaling to 0-1 on test data) of original data
minValues = min(allWeights,[],1);
maxValues = max(allWeights,[],1);
compFactor = 1./(maxValues-minValues);
save(pT.normValues, 'minValues', 'maxValues', 'compFactor');

% Normalize allWeights 
allWeights = bsxfun(@minus,allWeights,minValues);
allWeights = bsxfun(@times,allWeights,compFactor);

% skeleton annotation (assumes dense label in bboxSmall)
trainingData = allWeights(labelIdx,:);
trainingLabels = allLabels;
save(pT.initalGroundTruth, 'trainingData', 'trainingLabels');

end

