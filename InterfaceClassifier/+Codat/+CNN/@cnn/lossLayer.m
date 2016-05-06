function [ loss, err, delta ] = lossLayer( cnet, prediction, target, targetWeights )
%LOSSLAYER Calculate prediction loss and output layer deltas.
% INPUT prediction: 3d or 4d cube of voxel-predictions (activity{end} from
%                   forwardPass)
%       target: Target output cube of same size as prediction.
%       targetWeights: Array of single of same size as prediction or a
%           single number specifying a weight for each output pixel, i.e.
%           the loss at the pixel is multiplied with that weight. The loss
%           function is normalized by the sum over all weights. Use a zero
%           weight to "mask" a pixel.
%
% NOTE This function normalizes the loss by the sum of the target weights.
%
% NOTE TargetMask should be 0 and 1 only for sofmax error functions.
%      Weighted inputs are not implemented in this case. Furthermore, the
%      matrices targetWeights(:,:,:,i) should be equal for all i when using
%      softmax.
%
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

err = prediction - target;
if isscalar(targetWeights)
    sTw = targetWeights.*numel(prediction);
else
    sTw = sum(targetWeights(:));
end
err = err.*targetWeights./sTw;

switch cnet.lossFunction
    case 'squared'
        loss = 0.5*sum(targetWeights(:).*(prediction(:) - target(:)).^2);
        delta = cnet.nonLinearityD{end}(prediction).*err;
    case 'cross-entropy'
        % warning: only implemented for sigmoid activation in last layer
        %clipping prediction for loss calculation (does not affect the
        %gradient calculation which is based on error - should be replaced
        %by a more optimized version at some point e.g. sofplus)
        prediction(prediction < 1e-7) = 1e-7;
        prediction(prediction > (1 - 1e-7)) = 1 - 1e-7;
        loss = target.*log(prediction) + (1 - target).*log(1 - prediction);
        loss = loss.*targetWeights;
        loss = -sum(loss(:));
        delta = err;
    case 'softmax'
        %warning: only implemented for sofmax activation in the last layer
        %clipping prediction for loss calculation (does not affect the
        %gradient calculation which is based on error - should be replaced
        %by a more optimized version at some point e.g. softplus)
        prediction(prediction < 1e-7) = 1e-7;
        prediction(prediction > 1 - 1e-7) = 1 - 1e-7;
        loss = target.*log(prediction);
        loss = loss.*targetWeights;
        loss = -sum(loss(:));
        delta = err;
    otherwise
        error('Loss function not implemented');
end

% normalize by number of labels
loss = loss./sTw;

if cnet.l2WeightDecayLambda > 0
    l2Regularizer = 0;
    for lyr = 2:cnet.layer
        l2Regularizer = l2Regularizer + sum(cnet.W{lyr}(:).^2);
    end
    loss = loss + cnet.l2WeightDecayLambda/2*l2Regularizer;
end


end
