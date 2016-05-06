classdef denseConvLayerTest < matlab.unittest.TestCase
    %DENSECONVLAYERTEST Test functionaly of dense convolutional layer.
    
    properties
        epsilon = 1e-7;
        maxDiff = 1e-4;
    end
    
    methods (Test)
        function testConvAlgs(testCase)
            %without sparsity
            layer = Codat.Net.Layer.denseConvolutionalLayer([3, 3, 3], [1, 1, 1], 2, 3, [1, 1, 1]);
            layer = layer.initializeWeights('fan-in',0,'tanh');
            layer = layer.castParams(@double);
            layer2 = layer.setConvAlg('fft2');
            layer3 = layer.setConvAlg('convn');
            input = randn(10,10,10,2);
            [~,Y1] = layer.forward(input);
            [~,Y2] = layer2.forward(input);
            [~,Y3] = layer3.forward(input);
            diff12 = max(abs(Y1(:) - Y2(:)));
            diff23 = max(abs(Y2(:) - Y3(:)));
            diff13 = max(abs(Y1(:) - Y3(:)));
            verifyEqual(testCase, testCase.maxDiff > max([diff12,diff23,diff13]), true);
            
            %with sparsity
            layer = Codat.Net.Layer.denseConvolutionalLayer([3, 3, 3], [1, 1, 1], 2, 3,[2, 2, 2]);
            layer = layer.initializeWeights('fan-in',0,'tanh');
            layer = layer.castParams(@double);
            layer2 = layer.setConvAlg('fft2');
            layer3 = layer.setConvAlg('convn');
            input = randn(10,10,10,2);
            [~,Y1] = layer.forward(input);
            [~,Y2] = layer2.forward(input);
            [~,Y3] = layer3.forward(input);
            diff12 = max(abs(Y1(:) - Y2(:)));
            diff23 = max(abs(Y2(:) - Y3(:)));
            diff13 = max(abs(Y1(:) - Y3(:)));
            verifyEqual(testCase, testCase.maxDiff > max([diff12,diff23,diff13]), true);
        end
        
        function testSetConvAlgAndParameterSetting(testCase)
            filterSize = [2, 2, 2; 3, 3, 3; 3, 3, 1];
            sparsityFactor = [1, 1, 1;2, 2, 2];
            for i = 1:size(filterSize,1)
                for j = 1:size(sparsityFactor,1)
                    layer = Codat.Net.Layer.denseDeconvolutionalLayer(filterSize(i,:), [1, 1, 1], 2, 3, sparsityFactor(j,:));
                    w = randn(layer.numParams,1,'single');
                    layer = layer.vec2Param(w);
                    vec = layer.param2Vec();
                    verifyEqual(testCase, isequal(vec,w), true);
                    
                    layer = layer.setConvAlg('fft2');
                    vec = layer.param2Vec();
                    verifyEqual(testCase, isequal(vec,w), true);
                    
                    layer = layer.setConvAlg('fft');
                    vec = layer.param2Vec();
                    verifyEqual(testCase, isequal(vec,w), true);
                    
                    layer = layer.setConvAlg('convn');
                    vec = layer.param2Vec();
                    verifyEqual(testCase, isequal(vec,w), true);
                end
            end
        end
    end
    
end

