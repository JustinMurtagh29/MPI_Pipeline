clc;
cd C:\Users\mberning\Desktop\active;

%% load part of the data
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
        [weights{sub2ind([3 3],x,y)}, weightLabels] = featureDesign(rawT, affT, seg, edges{sub2ind([3 3],x,y)});
        time = toc;
        display([num2str(ceil(time/60)) ' min']);
    end
end
save('C:\Users\mberning\Desktop\active\wip.mat');

%% reformat & save (step(s) above will take >2h on the small traning dataset)
load('C:\Users\mberning\Desktop\active\wip.mat');
clear x y trueEdges time raw* aff* seg*;
edges = cell2mat(edges);
labelEdges = cell2mat(labelEdges);
weights = cell2mat(weights);
% add constant weight
weights(:,end+1) = ones(1,size(weights,1));
weightLabels{end+1} = 'constant feature';
% Split into training and test data and rename for convenience
testSize = 2000;
train.X = weights(1:end-testSize,:)';
train.y = labelEdges(1:end-testSize)';
train.edges = edges(1:end-testSize,:)';
test.X = weights(end-testSize+1:end,:)';
test.y = labelEdges(end-testSize+1:end);
test.edges = edges(end-testSize+1:end,:);
clear edges labelEdges testSize weights;
save('C:\Users\mberning\Desktop\active\beforeGaussianRegression.mat');

%% visualize features
visualizeEdgeFeatures(train, test, weightLabels);

%% with normalizing weights (without constant feature)
load('C:\Users\mberning\Desktop\active\beforeGaussianRegression.mat');
w.Mu = zeros(size(train.X,1),1);
w.Sigma = .3*eye(size(train.X,1));

meanX = mean([train.X(1:223,:), test.X(1:223,:)],2);
stdX = std([train.X(1:223,:), test.X(1:223,:)],[],2);
train.X(1:223,:) = train.X(1:223,:) - repmat(meanX,1,size(train.X,2));
train.X(1:223,:) = train.X(1:223,:) ./ repmat(stdX,1,size(train.X,2));
test.X(1:223,:) = test.X(1:223,:) - repmat(meanX,1,size(test.X,2));
test.X(1:223,:) = test.X(1:223,:) ./ repmat(stdX,1,size(test.X,2));
[p,ef, vf, approxLogMarginalLikelihood] = binaryGPclassifier(w, train.X, train.y', test.X');

display(['approximate marginal log likelihood: ' num2str(approxLogMarginalLikelihood)]);

%% plot prediction compared to ground truth on training set
figure;
[a, idx] = sort(p);
plot(a);
hold on;
b = test.y(idx);
d = sum(reshape(b, 100, 20),1)/100;
bar(50:100:1950,d);
xlabel('test case vector resorted by probabilty predicted by classifier for class 1');
legend('probabilty of belonging to class 1 predicted by classifier', 'frequency of class 1 label in bin according to ground truth', 'Location', 'Best');
title('Careful: Ground truth here is still comparison of manuell and automated segmentation (no skeletons yet)');

%% plot mean and variance of latent function
figure;
plot((1:2000)',ef(idx),'k','LineWidth',2);
hold on;
plot((1:2000)',ef(idx)-2*vf(idx),':k','LineWidth',2);
plot((1:2000)',ef(idx)+2*vf(idx),':k','LineWidth',2);

%% error reject tradeoff
figure;
p = p;
testP = .5:0.01:.99;
for i=1:length(testP)
    classified1 = p > testP(i);
    classified0 = p < (1 - testP(i));
    rejected = p > (1 - testP(i)) & p < testP(i);
    percentageRejected(i) = sum(rejected)/numel(rejected);
    percentageCorrectClassification(i) = (sum(classified1(classified1) == test.y(classified1)) + ...
        sum(~classified0(classified0) == test.y(classified0)))/sum(~rejected);
end
plot(100*percentageRejected, 100*percentageCorrectClassification);
xlabel('tasks rejected [%]');
ylabel('correctly classified [%]');
title('Error-Rejection tradeoff');

%% ROC 
figure;
p = p;
testP = 0:0.01:.99;
for i=1:length(testP)
    classified1 = p > testP(i);
    classified0 = p < testP(i);
    truePositiveRate(i) = sum(classified1(classified1) == test.y(classified1))/sum(test.y);
    falsePositiveRate(i) = sum(classified1(classified1) ~= test.y(classified1))/sum(~test.y);
end
plot(100*falsePositiveRate, 100*truePositiveRate);
xlabel('false positive rate [%]');
ylabel('true positive rate [%]');
title('ROC');

%% without normalizing weights
load('C:\Users\mberning\Desktop\active\beforeGaussianRegression.mat');
w.Mu = zeros(size(train.X,1),1);
w.Sigma = diag(var(train.X, [], 2));

[p, ef, vf, approxLogMarginalLikelihood] = binaryGPclassifier(w, train.X, train.y', test.X);

%% 
load('C:\Users\mberning\Desktop\active\beforeGaussianRegression.mat');
w.Mu = zeros(1);
w.Sigma = ones(1);

dimToTest = 17;
meanX = mean([train.X(dimToTest,:), test.X(dimToTest,:)],2);
stdX = std([train.X(dimToTest,:), test.X(dimToTest,:)],[],2);
train.X = train.X(17,:) - repmat(meanX,1,size(train.X,2));
train.X = train.X ./ repmat(stdX,1,size(train.X,2));
test.X = test.X(17,:) - repmat(meanX,1,size(test.X,2));
test.X = test.X ./ repmat(stdX,1,size(test.X,2));

figure;
for i=2:length(train.X)
    [p, ef, vf, approxLogMarginalLikelihood] = binaryGPclassifier(w, train.X(:,1:i), train.y(:,1:i)', test.X');
    [a,idx] = sort(test.X);
    plot(test.X(idx), p(idx), '-b');
    hold on;
    plot(train.X(:,1:i), train.y(:,1:i), 'xr');
    hold off;
    pause(1);
end