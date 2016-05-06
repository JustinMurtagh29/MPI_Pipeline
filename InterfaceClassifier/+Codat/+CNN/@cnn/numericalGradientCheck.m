function  [diff,numgrad,grad] = numericalGradientCheck( cnet )
%NUMERICALGRADIENTCHECK Perform a numerical gradient check by comparing the
% result of backprop to finite differences.
% INPUT options: (Optional) Struct containing the following field.
%           'testWeighted': Logical specifying whether the targetMask
%               should be logical or positive weights.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

rng('shuffle');
cnet.l2WeightDecayLambda = 1;
cnet.dropout = 0.5.*ones(cnet.layer,1);
cnet.isTraining = true;
cnet = cnet.setParamsToActvtClass(@double); %otherwise numerical gradient will most likely be incorrect
targetSize = [5 5 5 cnet.featureMaps(end)];
input = randn([cnet.border + targetSize(1:3) cnet.featureMaps(1)]);
target = randn(targetSize);
targetWeights = rand(size(target));
if strcmp(cnet.lossFunction,'cross-entropy')
    target(target < 0) = 0;
    target(target > 1) = 1;
elseif strcmp(cnet.lossFunction,'softmax')
    target = cnet.softmax(target);
    targetWeights = repmat((targetWeights(:,:,:,1) - 0.5) > 0,1,1,1,cnet.featureMaps(end));
end
epsilon = 1e-6;
j = 1;
rng(1);
[ activity,dropoutMask,mpInd, bn ] = forwardPass( cnet, input );
gradient = backprop( cnet,activity,dropoutMask,mpInd,bn,target,targetWeights );

for lyr = 2:cnet.layer
    %bias check
    for fm = 1:cnet.featureMaps(lyr)
        NN1 = cnet;
        NN2 = cnet;
        NN1.b{lyr}(fm) = NN1.b{lyr}(fm) + epsilon;
        NN2.b{lyr}(fm) = NN2.b{lyr}(fm) - epsilon;
        rng(1);
        activity = forwardPass( NN1, input );
        pos = NN1.lossLayer(activity{end}, target, targetWeights );
        rng(1);
        activity = forwardPass( NN2, input );
        neg = NN2.lossLayer(activity{end}, target, targetWeights );

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
            activity = forwardPass( NN1, input );
            pos = NN1.lossLayer(activity{end}, target, targetWeights );
            rng(1);
            activity = forwardPass( NN2, input );
            neg = NN2.lossLayer(activity{end}, target, targetWeights );
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
            activity = forwardPass( NN1, input );
            pos = NN1.lossLayer(activity{end}, target, targetWeights );
            rng(1);
            activity = forwardPass( NN2, input );
            neg = NN2.lossLayer(activity{end}, target, targetWeights );

            numgrad(j) = (pos-neg)/(2*epsilon);
            grad(j) = gradient.W{lyr}(idx);
            j = j + 1;
        end
    end
    
    %bn check
    if cnet.batchNorm(lyr)
        for i = 1:numel(cnet.bn_beta{lyr})
            NN1 = cnet;
            NN2 = cnet;
            NN1.bn_beta{lyr}(i) = NN1.bn_beta{lyr}(i) + epsilon;
            NN2.bn_beta{lyr}(i) = NN2.bn_beta{lyr}(i) - epsilon;
            rng(1);
            activity = forwardPass( NN1, input );
            pos = NN1.lossLayer(activity{end}, target, targetWeights );
            rng(1);
            activity = forwardPass( NN2, input );
            neg = NN2.lossLayer(activity{end}, target, targetWeights );
            numgrad(j) = (pos-neg)/(2*epsilon);
            grad(j) = gradient.beta{lyr}(i);
            j = j + 1;
        end
        
        for i = 1:numel(cnet.bn_gamma{lyr})
            NN1 = cnet;
            NN2 = cnet;
            NN1.bn_gamma{lyr}(i) = NN1.bn_gamma{lyr}(i) + epsilon;
            NN2.bn_gamma{lyr}(i) = NN2.bn_gamma{lyr}(i) - epsilon;
            rng(1);
            activity = forwardPass( NN1, input );
            pos = NN1.lossLayer(activity{end}, target, targetWeights );
            rng(1);
            activity = forwardPass( NN2, input );
            neg = NN2.lossLayer(activity{end}, target, targetWeights );
            numgrad(j) = (pos-neg)/(2*epsilon);
            grad(j) = gradient.gamma{lyr}(i);
            j = j + 1;
        end
    end
end

diff = max(abs(numgrad(:) - grad(:)));
tmp = [numgrad; grad];

end

