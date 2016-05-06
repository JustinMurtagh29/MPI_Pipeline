function [ activity, dropoutMask, mpInd, bn] = forwardFromTo( cnet, start_layer, end_layer, input )
%FORWARDFROMTO Forward propagation.
% INPUT layer_start: Index of layer to start forward pass.
%       layer_end: Index of target layer.
%       input: 4d cube where the first three dimensions are the data and
%              the forth dimension is the number of input channels/feature
%              maps in the corresponding layer.
% OUTPUT activity: Cell array index by layer where
%                  each entry contains the activity of the corresponding
%                  layer.
%        dropoutMask: Cell array indexed by layer of the dropout mask
%                     applied to the activities for backpropagation.
%        mpInd: Indices from max pooling.
%        bn: Cell array of length cnet.layer of quantities required for
%           batch normalization backpropagation. Each row corresponds to
%           one layer and contains the input to the bn layer, bn_beta,
%           bn_gamma, bn_muInf, bn_sig2Inf.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

%initialize cells
input = cnet.actvtClass(input);
activity = cell(cnet.layer,1);
dropoutMask = cell(1,cnet.layer);
bn = cell(cnet.layer,5);
mpInd = cell(1,cnet.layer);
activity{1} = input;

%dropout for input if specified
if start_layer == 1 && cnet.dropout(1) > 0 && cnet.isTraining
    if isa(input,'gpuArray')
        dM = gpuArray.rand(cnet.mSize(input,1:3)) < (1 - cnet.dropout(1));
    else
        dM = cnet.actvtClass(binornd(1,1 - cnet.dropout(1),cnet.mSize(input,1:3)));
    end
    activity{1} = bsxfun(@times,activity{1},dM)./(1-cnet.dropout(1));
end

%iterate over layers
for lyr = max(start_layer,2):min(cnet.layer, end_layer)
    %conv layer
    activityWithoutNL = cnet.convLayerFwd(activity{lyr - 1},lyr);

    %mp layer
    if cnet.maxPool(lyr)
        [activityWithoutNL,mpInd{lyr}] = cnet.maxPooling(activityWithoutNL,lyr);
    end

    %shortcut in
    if cnet.shortcut(lyr) > 0
        activityWithoutNL = activityWithoutNL + cnet.cropActivation(activity{cnet.shortcut(lyr)},size(activityWithoutNL));
    end
    
    %batch normalization
    if cnet.batchNorm(lyr) && lyr < cnet.layer
        bn{lyr,1} = activityWithoutNL;
        [activityWithoutNL, bn{lyr,2:end}] = Codat.NN.batchNormalization(activityWithoutNL, ...
            cnet.bn_beta{lyr}, cnet.bn_gamma{lyr}, cnet.bn_muInf{lyr}, cnet.bn_sig2Inf{lyr}, cnet.isTraining);
    end

    %dropout
    if cnet.dropout(lyr) > 0 && cnet.isTraining && lyr < cnet.layer
        if isa(activityWithoutNL,'gpuArray')
            dropoutMask{lyr} = gpuArray.rand(cnet.mSize(activityWithoutNL,1:3)) < (1 - cnet.dropout(lyr));
        else
            dropoutMask{lyr} = cnet.actvtClass(binornd(1,1 - cnet.dropout(lyr),cnet.mSize(activityWithoutNL,1:3)));
        end
        activityWithoutNL = bsxfun(@times,activityWithoutNL,dropoutMask{lyr})./(1-cnet.dropout(lyr));
    end

    %non-linearity
    activity{lyr} = cnet.nonLinearity{lyr}(activityWithoutNL);
end

end
