function [ gradient,delta ] = backwardFromTo( cnet, start_layer, end_layer, delta, activity, mpInd, dropoutMask, bn, applyOutputNonLinearity )
%BACKWARDFROMTO Backpropagation of errors through cnet.
% INPUT start: Starting layer for backpropagation. (Must be cnet.layer or
%              smaller)
%       end: End layer where error backpropagation is stopped. (Must be
%            1 or larger).
%       delta: Errors at start layer. Must have same layer as the
%           activity from forward propagation in the respective layer.
%       activity: Activity from forward pass.
%       mpInd: Max-pooling indices from forward propagation.
%       dropoutMask: Dropout mask used during forward pass. Set
%           cnet.isTraining to false if no should be used.
%       applyOutputNonLinearity: (Optional) Logical specifying whether
%           backpropgation is performed through the non-linearity in the
%           last layer (this is done in the loss function by default).
%           (Default: false)
%       bn: Batch normalization information (see forward pass).
% OUTPUT delta: Error at end layer.
%        gradient: Gradient for parameters of all layers.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if ~exist('applyOutputNonLinearity','var') || isempty(applyOutputNonLinearity)
    applyOutputNonLinearity = false;
end

start_layer = min(start_layer,cnet.layer);
end_layer = max(end_layer, 1);

%Initialize gradient
gradient.W = cell(cnet.layer,1);
gradient.b = cell(cnet.layer,1);
gradient.beta = cell(cnet.layer,1);
gradient.gamma = cell(cnet.layer,1);
shortcutDelta = cell(cnet.layer,1);

%backprop through hidden layers
prop_down = true;
for lyr = start_layer:-1:max(end_layer,2)
    %add delta from shortcut connection to higher layer
    if ~isempty(shortcutDelta{lyr}) && lyr > 1
        delta = delta + Codat.CNN.cnn.padArray(shortcutDelta{lyr},size(delta));
        shortcutDelta{lyr} = [];
    end

    %backprop through non-linearity
    if lyr < cnet.layer || applyOutputNonLinearity
        delta = delta.*cnet.nonLinearityD{lyr}(activity{lyr});
    end

    %backprop through dropout layer
    if cnet.dropout(lyr) > 0 && lyr < cnet.layer
        delta = bsxfun(@times,delta,dropoutMask{lyr})./(1-cnet.dropout(lyr));
    end
    
    %backprop through batch normalization layer
    if cnet.batchNorm(lyr) && lyr < cnet.layer
        [ delta, gradient.beta{lyr}, gradient.gamma{lyr} ] = Codat.NN.batchNormalizationBwd( delta, ...
            bn{lyr,1}, bn{lyr,2}, bn{lyr,3}, cnet.bn_gamma{lyr} );
    end

    %send delta via shortcut connection to lower layer
    if cnet.shortcut(lyr) > 0
        shortcutDelta{cnet.shortcut(lyr)} = delta;
    end

    %backprop through max pooling layer
    if cnet.maxPool(lyr)
        delta = cnet.maxPoolingBwd(delta,mpInd{lyr},lyr);
    end

    %backpropagate through convolution
    if lyr == 2 && end_layer == 2
        %dont propagate delta to first layer if not necessary
        prop_down = false;
    end
    [delta, gradient.W{lyr}, gradient.b{lyr}] = cnet.convLayerBwd( ...
        activity{lyr - 1}, delta, lyr, prop_down);
end

end
