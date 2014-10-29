% Load manuell segmentation
load /zdata/manuel/data/cortex/activeTraining/trainingSet.mat;

% Construct graph for training data & find ground truth
for x=1:3
    for y=1:3
        bbox = [2001+(x-1)*100 2000+x*100,...
		2001+(y-1)*100 2000+y*100,...
		2001:2100];
	temp = graphConstruction(bbox); 
	seg = loadSegData(temp.seg.root, temp.seg.prefix, temp.seg.bboxSmall);
       	load(temp.edgeFile);
	load(temp.weightFile);
	trueEdges = approximateGroundTruth(seg, seg_manuell(bbox));
	train.X{sub2ind([3 3],x,y)} = weights;
	train.y{sub2ind([3 3],x,y)} = findLabels(edges{sub2ind([3 3],x,y)}, trueEdges );
	train.edges{sub2ind([3 3],x,y)} = edges;
	trainingData{sub2ind([3 3],x,y)}.parameter = temp;
    end
end

% Concatenate results
train.X = cell2mat(train.X);
train.y = cell2mat(train.y);

% Load test data & put into structure
load(paramBG.edgeFile);
load(paramBG.weightFile);
test.X = weights;
test.edges = edges;

% Normalize weights to statistics of training data
nrFeature = size(train.X,2);
test.X = test.X - repmat(mean(train.X,2),1,nrFeature);
test.X = test.X ./ repmat(std(train.X,[],2),1,nrFeature);
train.X = train.X - repmat(mean(train.X,2),1,nrFeature);
train.X = train.X ./ repmat(std(train.X,[],2),1,nrFeature);

% add constant weight
test.X(:,end+1) = ones(1,size(test.X,1));
train.X(:,end+1) = ones(1,size(train.X,1));

% define prior and classify
w.Mu = zeros(size(train.X,1),1);
w.Sigma = .3*eye(size(train.X,1));
[p, latMean, latVar, nlml] = binaryGPclassifier(w, train.X, train.y', test.X);


