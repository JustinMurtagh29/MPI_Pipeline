clc;
cd /home/mberning/code/active/;
matlabpool 4;

%% load training data
profile on;
weights = cell(9,1);
edges = cell(9,1);
labelEdges = cell(9,1);
load dataWithGroundTruth.mat;
raw = single(raw);
raw = raw - mean(raw(:));
raw = raw ./ std(raw(:));
for x=1:3
    for y=1:3
        display([num2str(1+(x-1)*100) ' ' num2str(x*100) ' ' num2str(1+(y-1)*100) ' ' num2str(y*100)]);
        tic;
        affT = aff(1+(x-1)*100:x*100,1+(y-1)*100:y*100,1:100);
        rawT = raw(1+(x-1)*100:x*100,1+(y-1)*100:y*100,1:100);
        seg_manuellT = seg_manuell(1+(x-1)*100:x*100,1+(y-1)*100:y*100,1:100);
        clear seg;
        %%calculate edges, label for edges and weights
        seg = redoSegmentation( affT );
        trueEdges = approximateGroundTruth(seg, seg_manuellT);
        edges{sub2ind([3 3],x,y)} = findEdges(seg);
        labelEdges{sub2ind([3 3],x,y)} = findLabels( edges{sub2ind([3 3],x,y)}, trueEdges );
        [weights{sub2ind([3 3],x,y)}, ~] = featureDesign(rawT, affT, seg, edges{sub2ind([3 3],x,y)});
        time = toc;
        display([num2str(ceil(time/60)) ' min']);
    end
end

load forMoritz.mat;
aff = aff(1:500,1:500,1:500);
raw = single(raw(1:500,1:500,1:500));
raw = raw - mean(raw(:));
raw = raw ./ std(raw(:));
display('large test cube');
tic;
%%calculate edges, label for edges and weights
clear seg;
seg = redoSegmentation( aff );
edges2 = findEdges(seg);
[weights2, weightLabels] = featureDesign(raw, aff, seg, edges2);
time = toc;
display([num2str(ceil(time/60)) ' min']);
save('wip.mat');
profile viewer;
profile off;

%% reformat & save (step(s) above will take >2h on the small traning dataset)
load('wip.mat');
clear x y trueEdges time rawT affT segT;
edges = cell2mat(edges);
labelEdges = cell2mat(labelEdges);
weights = cell2mat(weights);
edges2 = cell2mat(edges2);
weights2 = cell2mat(weights2);
% add constant weight
weights(:,end+1) = ones(1,size(weights,1));
weights2(:,end+1) = ones(1,size(weights2,1));
weightLabels{end+1} = 'constant feature';
% Split into training and test data and rename for convenience
train.X = weights;
train.y = labelEdges;
train.edges = edges;
testSize = 2000;
train.X = weights(1:end-testSize,:)';
train.y = labelEdges(1:end-testSize)';
train.edges = edges(1:end-testSize,:)';
test.X = weights(end-testSize+1:end,:)';
test.y = labelEdges(end-testSize+1:end);
test.edges = edges(end-testSize+1:end,:);
clear edges labelEdges testSize weights;
save('beforeGaussianRegression.mat');

%% visualize features
visualizeEdgeFeatures(train, test, weightLabels);

%% with normalizing weights (without constant feature)
load('beforeGaussianRegression.mat');
w.Mu = zeros(size(train.X,1),1);
w.Sigma = .3*eye(size(train.X,1));

train.X(1:223,:) = train.X(1:223,:) - repmat(mean(train.X(1:223,:),2),1,size(train.X,2));
train.X(1:223,:) = train.X(1:223,:) ./ repmat(std(train.X(1:223,:),[],2),1,size(train.X,2));
test.X(1:223,:) = test.X(1:223,:) - repmat(mean(train.X(1:223,:),2),1,size(test.X,2));
test.X(1:223,:) = test.X(1:223,:) ./ repmat(std(train.X(1:223,:),[],2),1,size(test.X,2));
[predictiveProbabilities, predictiveMean, predictiveVariance, approxLogMarginalLikelihood] = binaryGPclassifier(w, train.X, train.y', test.X);

display(['Negative Log Likelihood: ' num2str(approxLogMarginalLikelihood)]);

[a, idx] = sort(predictiveProbabilities);
plot(a);
hold on;
b = test.y(idx);
d = sum(reshape(b, 100, 20),1)/100;
bar(50:100:1950,d);
xlabel('test case vector resorted by probabilty predicted by classifier for class 1');
legend('probabilty of belonging to class 1 predicted by classifier', 'frequency of class 1 label in bin according to ground truth', 'Location', 'Best');
title('Careful: Ground truth here is still comparison of manuell and automated segmentation (no skeletons yet)');

plot((1:2000)',predictiveMean(idx),'k','LineWidth',2);
hold on;
plot((1:2000)',predictiveMean(idx)-2*predictiveVariance(idx),':k','LineWidth',2);
plot((1:2000)',predictiveMean(idx)+2*predictiveVariance(idx),':k','LineWidth',2);
bar(50:100:1950,d);

%% without normalizing weights
load('beforeGaussianRegression.mat');
w.Mu = zeros(size(train.X,1),1);
w.Sigma = diag(var(train.X, [], 2));

[predictiveProbabilities, predictiveMean, predictiveVariance, approxLogMarginalLikelihood] = binaryGPclassifier(w, train.X, train.y', test.X);

