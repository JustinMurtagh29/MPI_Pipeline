function [ activity, dropoutMask, mpInd, bn ] = forwardPass( cnet, input )
%FORWARDPASS Forward pass through cnet.
% INPUT input: 4d cube where the first three dimensions are the data and
%              the forth dimension is the number of input channels.
% OUTPUT activity: Cell array index by layer where
%                  each entry contains the activity of the corresponding
%                  layer.
%        dropoutMask: Cell array indexed by layer of the dropout mask
%                     applied to the activities for backpropagation.
%        mpInd: Indices from max pooling.
%        bn: Cell array of length cnet.layer of quantities required for
%           batch normalization backpropagation. The cell for each layer
%           contains five further cells with the input to the bn during
%           forward pass, bn beta, bn gamma, bn_muInf, bn_sig2Inf.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

[activity, dropoutMask, mpInd, bn] = cnet.forwardFromTo(2, cnet.layer, input);

end
