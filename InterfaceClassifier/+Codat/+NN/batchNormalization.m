function [ X, muB, sig2B, muInf, sig2Inf ] = batchNormalization( X, beta, gamma, muInf, sig2Inf, isTraining, alpha, epsilon )
%BATCHNORMALIZATION Perform batch-normalization on input.
% Performs batch normalization on 4d arrays assuming that the fourth
% dimension corresponds to features.
% INPUT X: 4d numerical array where the first three dimensions correspond
%           to spatial dimensions and the fourth to feature maps.
%       beta: [1x1x1xN] numerical array, where N = size(X,4). Bias
%       	parameter for the normalize input.
%       gamma: [1x1x1xN] numerical array, where N = size(X,4). Scaling
%       	parameter for the normalize input.
%       muInf: (Optional) [1x1x1xN] numerical array where N = size(X,4) of
%           feature map means of the batch normalization layer for
%           prediction. If this input is specified then muInf is updated
%           via an exponential moving average for usage during validation.
%           For prediction muInf should be calculated from the whole
%           training set.
%       isTraining: (Optional) Bool specifying if network is training.
%           (Default: true)
%       sig2Inf: (Optional) as muInf but for the variances of each feature
%           map.
%       alpha: (Optional) Coefficient of weight decrease in the exponential
%           moving average.
%           (Default: 0.5 if muInf and sig2Inf are specified).
%       epsilon: (Optional) Conditioning paramter.
%           (Default: 1e-7)
% OUTPUT X: The batch-normalized input.
%        muB: [1x1x1xN] numerical array where N = size(X,4) of means for
%           each feature map of X.
%        sigmaB: [1x1x1xN] numerical array where N = size(X,4) of variances
%        	for each feature map of X.
%        muInf: The updated muInf. (Only if muInf was specified as input.)
%        sig2Inf: The updated sig2Inf. (Only if muInf was specified as
%        	(Only if muInf was specified as input.)
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if ~exist('isTraining','var') || isempty(isTraining)
    isTraining = true;
end
if ~exist('epsilon','var') || isempty(epsilon)
    epsilon = 1e-7;
end

%actual batch normalization
if isTraining
    %calculate batch statistics during training
    sX = size(X);
    m = prod(sX(1:3));
    muB = Codat.Lib.sumdims(X,1:3)./m;
    sig2B = Codat.Lib.sumdims(bsxfun(@minus,X,muB).^2,1:3)./m;
else
    %use inference statistics during prediction
    muB = muInf;
    sig2B = sig2Inf;
end
X = bsxfun(@minus,X,muB);
X = bsxfun(@rdivide,X,sqrt(sig2B + epsilon));
X = bsxfun(@times,X,gamma);
X = bsxfun(@plus,X,beta);

%update inference parameter via exponential moving average
if isTraining && (exist('muInf','var') && ~isempty(muInf)) && (exist('sig2Inf','var') && ~isempty(sig2Inf))
    if ~exist('alpha','var') || isempty(alpha)
        alpha = 0.999;
    end
    muInf = alpha.*muInf + (1-alpha).*muB;
    sig2Inf = alpha.*sig2Inf + (1-alpha).*sig2B;
end

end
