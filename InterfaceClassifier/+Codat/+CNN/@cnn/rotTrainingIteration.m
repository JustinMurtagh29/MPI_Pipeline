function [ cnet, loss, prediction, gradient ] = rotTrainingIteration( cnet, input, target, targetMask )
%ROTTRAININGITERATION Similar to trainingIterations but applies the same
%cnet to rotated versions of the input.
% INPUT input: Input cube
%       target: Target cube
% NOTE There are no checks here, so check for the right types/classes of
% parameters and data, correct combination of activation and error
% functions and correct input/output sizes before you call this function.
%
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

cnet.isTraining = true;
activity = cell(4,1);
dropoutMask = cell(4,1);
mpInd = cell(4,1);
bn = cell(4,1);
for rotIter = 1:4
    [ activity{rotIter}, dropoutMask{rotIter}, mpInd{rotIter}, bn{rotIter} ] = forwardPass( cnet, rot90(input,rotIter) );
end
prediction = zeros(size(activity{1}{end}),'like',activity{1}{end});
for rotIter = 1:4
    prediction = prediction + rot90(activity{rotIter}{end},4-rotIter)./4;
end
[loss,delta] = cnet.lossLayer(prediction, target, targetMask );
l2decay = cnet.l2WeightDecayLambda;
for rotIter = 1:4
    currGradient = cnet.backwardFromTo( cnet.layer, 2, ...
        rot90(delta,rotIter)./4, activity{rotIter}, mpInd{rotIter}, ...
        dropoutMask{rotIter}, bn{rotIter}, true );

    if rotIter == 1
        gradient = cnet.param2Vec(currGradient);
        cnet.l2WeightDecayLambda = 0; %do weight decay only once
    else
        gradient = gradient + cnet.param2Vec(currGradient);
    end
end
cnet.l2WeightDecayLambda = l2decay;

paramsOld = cnet.param2Vec();
[paramsNew,optimizer] = cnet.optimizer.optimize(paramsOld,gradient);
cnet = cnet.vec2Param(paramsNew);
cnet.optimizer = optimizer;
loss = gather(loss);

end
