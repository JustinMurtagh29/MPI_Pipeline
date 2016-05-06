function  [diff,numgrad,grad] = numericalGradientCheckRot( cnet, options )
%NUMERICALGRADIENTCHECK Perform a numerical gradient check by comparing the
% result of backprop to finite differences.
% INPUT options: (Optional) Struct containing the following field.
%           'testWeighted': Logical specifying whether the targetMask
%               should be logical or positive weights.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if ~exist('options','var') || isempty(options)
    options = struct;
end

rng('shuffle');
cnet.l2WeightDecayLambda = 0;
cnet.dropout = 0*ones(cnet.layer,1);
cnet.isTraining = true;
cnet = cnet.setParamsToActvtClass(@double); %otherwise numerical gradient will most likely be incorrect
targetSize = [20 20 20 cnet.featureMaps(end)];
input = randn([cnet.border + targetSize(1:3) cnet.featureMaps(1)]);
target = randn(targetSize);
if isfield(options,'testWeighted') && options.testWeighted
    targetMask = rand(size(target));
else
    targetMask = logical(repmat(binornd(1,0.1,targetSize(1:3)),1,1,1,cnet.featureMaps(end)));
end
if strcmp(cnet.lossFunction,'cross-entropy')
    target(target < 0) = 0;
    target(target > 1) = 1;
elseif strcmp(cnet.lossFunction,'softmax')
    target = cnet.softmax(target);
end
epsilon = 1e-8;
j = 1;
rng(1);
[~,~,~,gradient] = cnet.rotTrainingIteration(input,target,targetMask);
gradient = cnet.vec2Param(gradient);

for lyr = 2:cnet.layer
    %bias check
    for fm = 1:cnet.featureMaps(lyr)
        NN1 = cnet;
        NN2 = cnet;
        NN1.b{lyr}(fm) = NN1.b{lyr}(fm) + epsilon;
        NN2.b{lyr}(fm) = NN2.b{lyr}(fm) - epsilon;
        rng(1);
        [~,pos] = NN1.rotTrainingIteration(input,target,targetMask);
        rng(1);
        [~,neg] = NN2.rotTrainingIteration(input,target,targetMask);
        numgrad(j) = (pos-neg)/(2*epsilon);
        grad(j) = gradient.b{lyr}(fm);
        j = j + 1;
    end
    
    %weight check
    if any(strcmp({'fft1','fft2'},cnet.convAlg))
        for i = 1:numel(cnet.W{lyr})
            NN1 = cnet;
            NN2 = cnet;
            NN1.W{lyr}(i) = NN1.W{lyr}(i) + epsilon;
            NN2.W{lyr}(i) = NN2.W{lyr}(i) - epsilon;
            rng(1);
            [~,pos] = NN1.rotTrainingIteration(input,target,targetMask);
            rng(1);
            [~,neg] = NN2.rotTrainingIteration(input,target,targetMask);
            
            numgrad(j) = (pos-neg)/(2*epsilon);
            grad(j) = gradient.W{lyr}(i);
            j = j + 1;
        end
    else
        %indices needed to calculate gradient for sparse weight matrix
        indices = find(repmat(cnet.Wmask{lyr},[1, 1, 1, cnet.featureMaps(lyr - 1), cnet.featureMaps(lyr)]));
        for i = 1:length(indices)
            NN1 = cnet;
            NN2 = cnet;
            idx = indices(i);
            NN1.W{lyr}(idx) = NN1.W{lyr}(idx) + epsilon;
            NN2.W{lyr}(idx) = NN2.W{lyr}(idx) - epsilon;

            rng(1);
            [~,pos] = NN1.rotTrainingIteration(input,target,targetMask);
            rng(1);
            [~,neg] = NN2.rotTrainingIteration(input,target,targetMask);

            numgrad(j) = (pos-neg)/(2*epsilon);
            grad(j) = gradient.W{lyr}(idx);
            j = j + 1;
        end
    end
end

diff = max(abs(numgrad(:) - grad(:)));
tmp = [numgrad; grad];
end

