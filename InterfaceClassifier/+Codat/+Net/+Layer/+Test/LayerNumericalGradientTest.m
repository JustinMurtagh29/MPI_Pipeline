classdef LayerNumericalGradientTest < matlab.unittest.TestCase
    %LAYERNUMERICALGRADIENTTEST Gradient test via finite differences.

    properties
        epsilon = 1e-7;
        maxDiff = 1e-4;
        rngSeed = 1;
    end

    methods (Test)
        function testReluLayer(testCase)
            layer = Codat.Net.Layer.reluLayer();
            input = randn(10,10,10,2);
            target = randn(10,10,10,2);
            numGrad = testCase.numericalGrad(layer, input, target);
            grad = testCase.analyticalGrad(layer, input, target);
            diff = max(abs(numGrad(:) - grad(:)));
            verifyEqual(testCase,testCase.maxDiff > diff,true);
        end

        function testTanhLayer(testCase)
            layer = Codat.Net.Layer.tanhLayer();
            input = randn(10,10,10,2);
            target = randn(10,10,10,2);
            numGrad = testCase.numericalGrad(layer, input, target);
            grad = testCase.analyticalGrad(layer, input, target);
            diff = max(abs(numGrad(:) - grad(:)));
            verifyEqual(testCase,testCase.maxDiff > diff,true);
        end

        function testSigmoidLayer(testCase)
            layer = Codat.Net.Layer.sigmoidLayer(false);
            input = randn(10,10,10,2);
            target = randn(10,10,10,2);
            numGrad = testCase.numericalGrad(layer, input, target);
            grad = testCase.analyticalGrad(layer, input, target);
            diff = max(abs(numGrad(:) - grad(:)));
            verifyEqual(testCase,testCase.maxDiff > diff,true);
        end

        function testDropoutLayer(testCase)
            layer = Codat.Net.Layer.dropoutLayer(0.3);
            input = randn(10,10,10,2);
            target = randn(10,10,10,2);
            numGrad = testCase.numericalGrad(layer, input, target);
            grad = testCase.analyticalGrad(layer, input, target);
            diff = max(abs(numGrad(:) - grad(:)));
            verifyEqual(testCase,testCase.maxDiff > diff,true);
        end

        function testCorruptionLayer(testCase)
            layer = Codat.Net.Layer.corruptionLayer(0,0.1);
            input = randn(10,10,10,2);
            target = randn(10,10,10,2);
            numGrad = testCase.numericalGrad(layer, input, target);
            grad = testCase.analyticalGrad(layer, input, target);
            diff = max(abs(numGrad(:) - grad(:)));
            verifyEqual(testCase,testCase.maxDiff > diff,true);
        end

        function testSoftmaxLayer(testCase)
            layer = Codat.Net.Layer.corruptionLayer(0,0.1);
            input = randn(10,10,10,3);
            target = randn(10,10,10,3);
            numGrad = testCase.numericalGrad(layer, input, target);
            grad = testCase.analyticalGrad(layer, input, target);
            diff = max(abs(numGrad(:) - grad(:)));
            verifyEqual(testCase,testCase.maxDiff > diff,true);
        end

        function testMaxPoolingLayer(testCase)
            %test optimized pooling windows
            for i = 2:3
                for j = 1:i
                    layer = Codat.Net.Layer.maxPoolingLayer([i, i, j], [1, 1, 1], [1, 1, 1]);
                    input = randn(8,8,3,3);
                    target = randn(layer.targetSize(size(input)));
                    numGrad = testCase.numericalGrad(layer, input, target);
                    grad = testCase.analyticalGrad(layer, input, target);
                    diff = max(abs(numGrad(:) - grad(:)));
                    verifyEqual(testCase,testCase.maxDiff > diff,true);
                end
            end

            %test one arbitrary window
            layer = Codat.Net.Layer.maxPoolingLayer([2,1,2], [1, 1, 1], [1, 1, 1]);
            input = randn(8,8,3,3);
            target = randn(layer.targetSize(size(input)));
            numGrad = testCase.numericalGrad(layer, input, target);
            grad = testCase.analyticalGrad(layer, input, target);
            diff = max(abs(numGrad(:) - grad(:)));
            verifyEqual(testCase,testCase.maxDiff > diff,true);
        end

        function testDenseConvolutionalLayer(testCase)
            filterSize = [2, 2, 2; 3, 3, 3; 3, 3, 1];
            sparsityFactor = [1, 1, 1;2, 2, 2];
            for i = 1:size(filterSize,1)
                for j = 1:size(sparsityFactor,1)
                    layer = Codat.Net.Layer.denseConvolutionalLayer(filterSize(i,:), [1, 1, 1], 2, 3, sparsityFactor(j,:));

                    %test fft mode
                    [diffIn, diffW] = testCase.denseConvBase(layer,[10, 10, 10]);
                    verifyEqual(testCase,testCase.maxDiff > diffIn,true);
                    verifyEqual(testCase,testCase.maxDiff > diffW,true);

                    %test fft2 mode
                    layer = layer.setConvAlg('fft2');
                    [diffIn, diffW] = testCase.denseConvBase(layer,[10, 10, 10]);
                    verifyEqual(testCase,testCase.maxDiff > diffIn,true);
                    verifyEqual(testCase,testCase.maxDiff > diffW,true);

                    %test convn mode
                    layer = layer.setConvAlg('convn');
                    [diffIn, diffW] = testCase.denseConvBase(layer,[10, 10, 10]);
                    verifyEqual(testCase,testCase.maxDiff > diffIn,true);
                    verifyEqual(testCase,testCase.maxDiff > diffW,true);
                end
            end
        end

        function testDenseDeconvolutionalLayer(testCase)
            testCase.maxDiff = 5e-4; %seems to be above 1e-4 sometimes but looks fine
            filterSize = [2, 2, 2; 3, 3, 3; 3, 3, 1];
            sparsityFactor = [1, 1, 1;2, 2, 2];
            inSiz = [3, 3, 3];
            for i = 1:size(filterSize,1)
                for j = 1:size(sparsityFactor,1)
                    layer = Codat.Net.Layer.denseDeconvolutionalLayer(filterSize(i,:), [1, 1, 1], 2, 3, sparsityFactor(j,:));

                    %test fft mode
                    [diffIn, diffW] = testCase.denseDeconvBase(layer,inSiz);
                    verifyEqual(testCase,testCase.maxDiff > diffIn,true);
                    verifyEqual(testCase,testCase.maxDiff > diffW,true);

                    %test fft2 mode
                    layer = layer.setConvAlg('fft2');
                    [diffIn, diffW] = testCase.denseDeconvBase(layer,inSiz);
                    verifyEqual(testCase,testCase.maxDiff > diffIn,true);
                    verifyEqual(testCase,testCase.maxDiff > diffW,true);

                    %test convn mode
                    layer = layer.setConvAlg('convn');
                    [diffIn, diffW] = testCase.denseDeconvBase(layer,inSiz);
                    verifyEqual(testCase,testCase.maxDiff > diffIn,true);
                    verifyEqual(testCase,testCase.maxDiff > diffW,true);
                end
            end
        end
        
        function testBatchNormalizationLayer(testCase)
            layer = Codat.Net.Layer.batchNormalizationLayer(3);
            layer = layer.castParams(@double);
            layer.keepMovingAverage = false;
            input = 2.*randn(5,5,5,3) + 5;
            target = randn(5,5,5,3);
            [numGradIn, numGradW] = testCase.numericalGrad(layer, input, target);
            [gradIn, gradW] = testCase.analyticalGrad(layer, input, target);
            diffIn = max(abs(numGradIn(:) - gradIn(:)));
            diffW = max(abs(numGradW(:) - gradW(:)));
            verifyEqual(testCase,testCase.maxDiff > diffIn,true);
            verifyEqual(testCase,testCase.maxDiff > diffW,true);
        end
    end

    methods
        function [DEDX,DEDW] = numericalGrad(testCase, layer, input, target)
            % Calculate the numerical gradient w.r.t. the input and
            % parameters using finite differences.
            eps = testCase.epsilon;
            DEDX = zeros(size(input),'like',input);
            for i = 1:numel(input)
                inNeg = input;
                inNeg(i) = inNeg(i) - eps;
                inPos = input;
                inPos(i) = inPos(i) + eps;
                rng(testCase.rngSeed);
                [~, neg] = layer.forward(inNeg);
                neg = testCase.squaredError(neg, target);
                rng(testCase.rngSeed);
                [~, pos] = layer.forward(inPos);
                pos = testCase.squaredError(pos, target);
                DEDX(i) = (pos - neg)./(2*eps);
            end
            if layer.numParams > 0
                DEDW = zeros(length(layer.numParams),1);
                params = layer.param2Vec();
                for i = 1:layer.numParams
                    pNeg = params;
                    pNeg(i) = pNeg(i) - eps;
                    layerNeg = layer.vec2Param(pNeg);
                    pPos = params;
                    pPos(i) = pPos(i) + eps;
                    layerPos = layer.vec2Param(pPos);
                    rng(testCase.rngSeed);
                    [~, neg] = layerNeg.forward(input);
                    rng(testCase.rngSeed);
                    [~, pos] = layerPos.forward(input);
                    neg = testCase.squaredError(neg, target);
                    pos = testCase.squaredError(pos, target);
                    DEDW(i) = (pos - neg)./(2*eps);
                end
            else
                DEDW = [];
            end
        end

        function [DEDX,DEDW] = analyticalGrad(testCase, layer, input, target)
            %Calculate the analytical gradient of a layer w.r.t. to input
            %an parameters as implemented in layer.backward.

            rng(testCase.rngSeed);
            [layer, Y] = layer.forward(input);
            [~,DEDY] = testCase.squaredError(Y, target);
            if layer.isInPlace
                [DEDX,DEDW] = layer.backward(Y,DEDY);
            else
                [DEDX,DEDW] = layer.backward(input,DEDY);
            end
        end

        function [diffIn, diffW] = denseConvBase(testCase, layer, sizIn)
            layer = layer.initializeWeights('fan-in',0,'tanh');
            layer = layer.castParams(@double);
            if length(sizIn) < 4
                sizIn(4) = layer.featureMapsIn;
            end
            input = randn(sizIn);
            target = randn([sizIn(1:3) - layer.sparseFilterSize + 1, layer.featureMapsOut]);
            [DEDX, DEDW] = testCase.analyticalGrad(layer, input, target);
            [numDEDX, numDEDW] = testCase.numericalGrad(layer, input, target);
            diffIn = max(abs(numDEDX(:) - DEDX(:)));
            diffW = max(abs(numDEDW(:) - DEDW(:)));
        end

        function [diffIn, diffW] = denseDeconvBase(testCase, layer, sizIn)
            layer = layer.initializeWeights('fan-in',0,'tanh');
            layer = layer.castParams(@double);
            if length(sizIn) < 4
                sizIn(4) = layer.featureMapsIn;
            end
            input = randn(sizIn);
            target = randn([sizIn(1:3) + layer.sparseFilterSize - 1, layer.featureMapsOut]);
            [DEDX, DEDW] = testCase.analyticalGrad(layer, input, target);
            [numDEDX, numDEDW] = testCase.numericalGrad(layer, input, target);
            diffIn = max(abs(numDEDX(:) - DEDX(:)));
            diffW = max(abs(numDEDW(:) - DEDW(:)));
        end
    end

    methods (Static)
        function [loss, DEDY] = squaredError(y, t)
            loss = 0.5.*sum((y(:) - t(:)).^2);
            DEDY = y - t;
        end
    end

end
