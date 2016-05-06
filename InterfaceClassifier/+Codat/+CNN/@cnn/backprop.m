function [ gradient, loss ] = backprop( cnet, activity, dropoutMask, mpInd, bn, target, targetWeights )
%BACKPROP Gradient calculation via error backpropagation.
% INPUT activity: Cell array of activities for each layer (first output of
%                 forward pass).
%       dropoutMask: Cell array of dropout masks used during forward
%                    propagation (second output of forward pass).
%       mpInd: Indices from max-pooling layers (third output of forward
%              pass).
%       target: Target output cube of same size as activity{end}.
%       targetMask: Mask for target voxels to discard. Use
%                   False(size(target)) to calculate the error for the
%                   whole target cube.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

%loss layer
[loss,~,delta] = cnet.lossLayer(activity{end}, target, targetWeights );

%parameter gradient
gradient = cnet.backwardFromTo( cnet.layer, 2, delta, activity, mpInd, ...
    dropoutMask, bn, false );

end
