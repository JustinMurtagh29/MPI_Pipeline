function [ DEDX, DEDb, DEDg ] = batchNormalizationBwd( DEDY, X, muB, sig2B, gamma, epsilon )
%BATCHNORMALIZATIONBWD Backward pass through batch-normalization layer.
% INPUT DEDY: 4d array containing the error derivative of the loss function
%           wrt to the output of the forward pass of the batch
%           normalization layer.
%       X: 4d array containing the input of the batch normalization layer
%           during forward pass.
%       muB: [1x1x1xN] numerical array containing muB output from forward
%           pass.
%       sig2B: [1x1x1xN] numerical array containing sig2B output from
%           forward pass.
%       gamma: [1x1x1xN] numerical array containing the current gamma.
%       epsilon: (Optional) 1d conditioning parameter.
%           (Default: 1e-7)
% OUTPUT DEDX: 4d array of backpropagated errors through.
%        DEDb: [1x1x1xN] numerical array  containing the gradient wrt to
%           the layer beta parameter.
%        DEDg: [1x1x1xN] numerical array containing the gradient wrt to
%           the layer gamma parameter.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if ~exist('epsilon','var') || isempty(epsilon)
    epsilon = 1e-7;
end

sX = size(X);
m = prod(sX(1:3));
X = bsxfun(@minus,X,muB);
DEDXh = bsxfun(@times,DEDY,gamma);
DEDsigB_2 = -0.5.*Codat.Lib.sumdims(DEDXh.*X,1:3).*(sig2B + epsilon).^(-3/2);
DEDmuB = -Codat.Lib.sumdims(DEDXh,1:3)./(sqrt(sig2B + epsilon)) - 2.*DEDsigB_2./m.*Codat.Lib.sumdims(X,1:3);
DEDX = bsxfun(@rdivide,DEDXh,sqrt(sig2B + epsilon)) + bsxfun(@plus,2.*bsxfun(@times,X,DEDsigB_2)./m,DEDmuB./m);
DEDg = Codat.Lib.sumdims(DEDY.*bsxfun(@rdivide,X,sqrt(sig2B + epsilon)),1:3);
DEDb = Codat.Lib.sumdims(DEDY,1:3);

end

