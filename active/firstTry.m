% Load manuell segmentation
%load /zdata/manuel/data/cortex/activeTraining/trainingSeg.mat;
% Load parameter from big graph construction
clc
load /zdata/manuel/data/cortex/activeTraining/testParam.mat;

% %% Construct graph for training data & find ground truth
% for x=1:3
%     for y=1:3
%         bbox = [2001+(x-1)*100 2000+x*100;...
% 		2001+(y-1)*100 2000+y*100;...
% 		2001 2100];
% 	temp = graphConstruction(bbox); 
% 	seg = loadSegData(temp.seg.root, temp.seg.prefix, temp.bboxSmall);
%     load(temp.edgeFile);
% 	load(temp.weightFile);
% 	trueEdges = approximateGroundTruth(seg, seg_manuell(bbox));
% 	train.X{sub2ind([3 3],x,y)} = real(weights);
% 	train.y{sub2ind([3 3],x,y)} = findLabels(edges{sub2ind([3 3],x,y)}, trueEdges );
% 	train.edges{sub2ind([3 3],x,y)} = edges;
% 	trainingData{sub2ind([3 3],x,y)}.parameter = temp;
%     end
% end
load /zdata/manuel/data/cortex/activeTraining/resultTrainingP1-377.mat

%% Concatenate results
train.X = real(weights');
train.y = labelEdges;

% Load test data & put into structure
load(paramBG.edgeFile);
load(paramBG.weightFile);
test.X = real(weights');
test.edges = edges;

% Normalize weights to statistics of training data
nrFeature = size(train.X,1);
test.X = test.X - repmat(mean(train.X,2),1,size(test.X,2));
test.X = test.X ./ repmat(std(train.X,[],2),1,size(test.X,2));
train.X = train.X - repmat(mean(train.X,2),1,size(train.X,2));
train.X = train.X ./ repmat(std(train.X,[],2),1,size(train.X,2));

% add constant weight
test.X(end+1,:) = ones(size(test.X,2),1);
train.X(end+1,:) = ones(size(train.X,2),1);

% define prior and classify
w.Mu = zeros(size(train.X,1),1);
w.Sigma = .3*eye(size(train.X,1));
[p, latMean, latVar, nlml] = binaryGPclassifier(w, train.X, train.y, test.X);

%% plot some stuff
figureLocation = '/zdata/manuel/sync/activeTraining/';
seg = loadSegData(paramBG.seg.root, paramBG.seg.prefix, paramBG.bboxSmall);
raw = loadRawData(paramBG.raw.root, paramBG.raw.prefix, paramBG.bboxSmall, 0);
weightLabel = generateWeightLabels();
%visualizeEdgeFeatures(train, test, weightLabel, figureLocation);
save([figureLocation 'beforeGP.mat'], 'raw', 'seg');

% make decision on error/reject of edges
errorCost = 1000; % cost of error in relation to rejection
[segNew, edgesNew, pNew] = joinSegments(seg, edges, p, errorCost);
save([figureLocation 'afterGP.mat'], 'raw', 'segNew');

% Resort to be more likely to look at connected pairs
[pResort, idx] = sort(pNew, 'descend');
edgesResort = edgesNew(idx,:);

% transfer everything to knowledgeDB
seg = loadSegData(paramBG.seg.root, paramBG.seg.prefix, paramBG.bboxBig);
seg(257:end-256,257:end-256,129:end-128) = segNew;
[settings, missions] = writeKnowledgeDB(paramBG, seg, edgesResort, pResort);

% visualization does not work with OpenGL harware acceleration (and thus is very slow)
%visualizeTasks(settings, missions, segNew);

