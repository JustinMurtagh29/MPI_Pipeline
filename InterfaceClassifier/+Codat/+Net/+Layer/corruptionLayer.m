classdef corruptionLayer < Codat.Net.Layer.baseLayer
    %CORRUPTIONLAYER Pointwise corruption with gaussian white noise.
    % Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>
    
    properties
        mean					%noise variance
        sigma					%noise standard deviation
    end
    
    methods
        function obj = corruptionLayer(mean, sigma)
            obj.mean = mean;
            obj. sigma = sigma;
            obj.isInPlace = true;
        end
        
        function [ layer, X ] = forward(layer, X)
            X = X + randn(size(X),'like',X);
        end
        
        function [DEDX, DEDW] = backward(~, ~, DEDY)
            DEDX = DEDY;
            DEDW = [];
        end
    end
    
end
