%% visualize features
load('C:\Users\mberning\Desktop\active\beforeGaussianRegression.mat');

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

%% calculate posterior distribution for w & predictions f* for regression
w2 = gaussianRegression( w, noiseSigmaSquared, X, Y);
plotMarginal2D(w2, [17 224], weightLabels);
f_star = makePredictions(w2, 17, x_star);
plotMarginal(f_star, 17, weightLabels);
f_star = makePredictions(w2, 224, x_star);
plotMarginal(f_star, 224, weightLabels);

%% switch to linear logistic classification

