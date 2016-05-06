classdef tanhLayer < Codat.Net.Layer.baseLayer
    %TANHLAYER Pointwise tanh non-linearity.
    % Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>
    
    properties
    end
    
    methods
        function obj = tanhLayer()
            obj.isInPlace = true;
        end
        
        function [ layer, X ] = forward(layer, X)
            X = 1.7159*tanh(0.66*X);
        end
        
        function [DEDX, DEDW] = backward(~, Y, DEDY)
            DEDX = DEDY.*0.66.*(1.7159 - 1./1.7159.*Y.^2);
            DEDW = [];
        end
    end
    
end

