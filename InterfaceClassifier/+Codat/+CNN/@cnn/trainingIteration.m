function [ cnet, loss, prediction ] = trainingIteration( cnet, input, target, targetWeights )
%TRAININGITERATION One training iteration of the cnn.
% INPUT input: Input cube
%       target: Target cube
% NOTE There are no checks here, so check for the right types/classes of
% parameters and data, correct combination of activation and error
% functions and correct input/output sizes before you call this function.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

cnet.isTraining = true;
[ activity, dropoutMask, mpInd, bn ] = forwardPass( cnet, input );

%set moving bn inference parameters to exp moving averages from training
cnet.bn_muInf = bn(:,4);
cnet.bn_sig2Inf = bn(:,5);

[ gradient, loss ] = backprop( cnet, activity, dropoutMask, mpInd, bn, target, targetWeights );
cnet = cnet.runOptimizer(gradient);
loss = gather(loss);
prediction = activity{end};

end

