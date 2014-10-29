clc;
cd C:\Users\mberning\Desktop\active;

%% load part of the data
load dataWithGroundTruth.mat;
aff = aff(1:100,1:100,1:100);
raw = single(raw);
raw = raw - mean(raw(:));
raw = raw ./ std(raw(:));
raw = raw(1:100,1:100,1:100);
seg_manuell = seg_manuell(1:100,1:100,1:100);
clear seg;

%% calculate edges, label for edges and weights
[seg, trueEdges] = approximateGroundTruth(aff, seg_manuell);
edges = findEdges(seg);
labelEdges = findLabels( edges, trueEdges );
[weights, weightLabels] = featureDesign(raw, aff, seg, edges);
% normalize weights
% weights = weights - repmat(mean(weights,1),size(weights,1),1);
% weights = weights ./ repmat(std(weights,0,1),size(weights,1),1);
% add constant weight
weights(:,end+1) = ones(1,size(weights,1));
weightLabels{end+1} = 'constant feature';
% save (step(s) above will take ~30min on the small dataset
save('C:\Users\mberning\Desktop\active\beforeGaussianRegression.mat');

%% visualize features
load('C:\Users\mberning\Desktop\active\beforeGaussianRegression.mat');
visualizeEdgeFeatures(weights, labelEdges, weightLabels);

%% define prior over weights and gaussian noise for observations
w.Mu = zeros(size(weights,2),1);
w.Sigma = .3*eye(size(weights,2));

noiseSigmaSquared = 0.2;
plotMarginal2D(w, [17 224], weightLabels);

%% set up variables for bayesian linear regression
X = weights';
Y = labelEdges';
x_star = -3:0.01:3;

%% calculate posterior distribution for w & predictions f*
w2 = gaussianRegression( w, noiseSigmaSquared, X, Y);
plotMarginal2D(w2, [17 224], weightLabels);
f_star = makePredictions(w2, 17, x_star);
plotMarginal(f_star, 17, weightLabels);
f_star = makePredictions(w2, 224, x_star);
plotMarginal(f_star, 224, weightLabels);

