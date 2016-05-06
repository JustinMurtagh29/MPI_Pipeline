classdef sigmoidLayer < Codat.Net.Layer.baseLayer
    %SIGMOIDLAYER Pointwise tanh non-linearity.
    %
    % PROPERTIES
    % onlyBackpass: Bool specifying whether error should be propagated back
    %       through non-linarity or simply passed back to the next layer. This
    %       options can be used if the softmax layer is followed by a softmax
    %       loss layer to simplify computations.
    %
    % NOTE This layer is mainly intended to be used directly before a
    % log loss layer with onlyBackpass set to true. For hidden layers it is
    % typically better to use a tanh layer or a relu layer.
    %     
    % Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>
    
    properties
        onlyBackpass = false;
    end
    
    methods
        function obj = sigmoidLayer(onlyBackpass)
            if exist('onlyBackpass','var') || ~isempty(onlyBackpass)
                obj.onlyBackpass = onlyBackpass;
            end
            obj.isInPlace = true;
        end
        
        function [ layer, X ] = forward(layer, X)
            X = sigmf(X,[1, 0]);
        end
        
        function [DEDX, DEDW] = backward(layer, Y, DEDY)
            if layer.onlyBackpass
                DEDX = DEDY;
            else
                DEDX = DEDY.*Y.*(1 - Y);
            end
            DEDW = [];
        end
    end
    
end
