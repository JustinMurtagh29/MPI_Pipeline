%% visualize features
load('C:\Users\mberning\Desktop\active\forPhilipp.mat');
visualizeEdgeFeatures(weights, labelEdges, weightLabels);

%% define prior over weights and constant gaussian noise for observations
w.Mu = zeros(size(weights,2),1);
w.Sigma = .3*eye(size(weights,2));
noiseSigmaSquared = 0.2;
% plot marginal distribution of w wrt 17th and 224th feature (17th
% classifciation result within border between objects, 224th constant feature (=1))
plotMarginal2D(w, [17 224], weightLabels);

%% set up variables for easier notation
X = weights';
Y = labelEdges';
x_star = -3:0.01:3;

%% calculate posterior distribution for w (regression)
w2 = gaussianRegression( w, noiseSigmaSquared, X, Y);
imagesc(w2.Sigma);
plotMarginal2D(w2, [17 224], weightLabels);
plotMarginal2D(w2, [91 96], weightLabels);
plotMarginal2D(w2, [51 56], weightLabels);

%% predictions f* (regression)
f_star = makePredictions(w2, 17, x_star);
plotMarginal(f_star, x_star, X(17,:), Y, weightLabels{17});
f_star = makePredictions(w2, 224, x_star);
plotMarginal(f_star, x_star, X(224,:), weightLabels{224});

%% different option
makePredictions_v2(w, X, Y, x_star, noiseSigmaSquared, 17);

%% switch to linear probit classification
w2 = gaussianProbit( w, noiseSigmaSquared, X', Y');
