%% load training data
load('C:\Users\mberning\Desktop\active\forPhilipp.mat');

%% set up variables for easier notation
X = weights';
Y = labelEdges;

%% define prior
w.Mu = zeros(size(weights,2),1);
w.Sigma = .3*eye(size(weights,2));

%% laplace approximation for binary bayesian classification
x_star = zeros(size(w.Mu,1),601);
x_star(1,:) = -3:0.01:3;
[p, mu, var, nlz] = binaryGPclassifier(w, X, Y, x_star);
