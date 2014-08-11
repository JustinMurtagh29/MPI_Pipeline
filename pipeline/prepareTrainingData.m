function prepareTrainingData(pT,mode)
% Collects normalized training data to pT.initalGroundTruth and values needed for normalization of test data to pT.normValues

allIdx = [];
allLabels = [];
allWeights = [];

if strcmp(mode,'edges')
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
    allWeights = normalizeDataForGP(allWeights, true, pT.gp.normValues);

elseif strcmp(mode,'glia')  
    for i=1:length(pT.local)
        % Load data of current densly skeletonized training region
        [~, seg] = loadSegData(pT.local(i).segFile, pT.tileBorder);
        load(pT.local(i).segmentFile);
        load(pT.local(i).segmentWeightFile);
        skel.file = pT.local(i).trainFile;
        skel.bbox = pT.local(i).bboxSmall;
        % Extract labels for all segments that intersect with any of the dense skeletons
        [labelIdx, labels] = extractGroundTruthFromNml_glia(seg, segments, skel);
        % Keep allWeights and just the labelled ones for normalization of test data & GP training respectively
        allWeights = [allWeights; segmentWeights];
        allIdx = [allIdx; labelIdx];
        allLabels = [allLabels; labels];
        clear seg segments segmentWeights skel labelIdx labels;
    end
    allWeights = normalizeDataForGP(allWeights, true, pT.glia.normValues);    
end


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

if strcmp(mode,'edges')
    save(pT.gp.initalGroundTruth, 'trainingData', 'trainingLabels', 'testData', 'testLabels');
    visualizeEdgeFeaturesNewest(pT);
elseif strcmp(mode,'glia')  
    save(pT.glia.initalGroundTruth, 'trainingData', 'trainingLabels', 'testData', 'testLabels');
end
end

