function X = augmentFeatures(X, dev)
%AUGMENTFEATURES Add random noise to features based on the feature mean.
% INPUT X: [NxM] float
%           Feature matrix. Rows correspond to instances and columns to
%           features.
%       dev: (Optional) float
%           Fraction of the feature mean that is used as the standard
%           deviation for random noise for the corresponding feature. i.e.
%           for each column the following is done
%           X(:,i) = X(:,i) + randn(length(X(:,i), 1), 1).*mean(X(:,i))*dev
%           (Default: 0.01)
% OUTPUT X: [NxM] float
%           The input features with added random noise.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if ~exist('dev', 'var') || isempty(dev)
    dev = 0.01;
end

m = mean(X, 1);
X = bsxfun(@plus, X, bsxfun(@times, randn(size(X), 'like', X), m.*dev));

end
