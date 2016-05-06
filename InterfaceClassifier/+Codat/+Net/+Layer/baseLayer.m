classdef baseLayer
    %LAYER Layer base class. Each layer should be derived from this class.
    %
    % PROPERTIES
    % numParams: Specifies the total number of learnable parameters of the
    %            layer.
    % stride: Vector of integer specifying the stride of the layer in each
    %         dimension. For dense pixel classification this stride is used
    %         to adapt the input window sizes of succeeding layer.
    % isInPlace: Bool specifying whether the layer performs its operations
    %            in place. If it performs the operation in place, then
    %            layer.backward should use the output Y of layer.forward in
    %            backpropagation. Otherwise it should use the input X of
    %            layer.forward.
    % Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

    properties
        %properties should be redefined appropriately
        numParams = 0
        stride = 1
        isInPlace
    end

    methods (Abstract)
        %forward: forward pass through layer
        [layer, Y] = forward(layer, X)

        %backward: backward pass through layer
        %The second input can be X or Y depending on whether to layer
        %performs in-place operations.
        [DEDX, DEDW] = backward(layer, XorY, DEDY)
    end

    methods
        %These methods should be adapted by each class containing
        %parameters. The implementation below is the standard for layers
        %without parameters.
        
        function layer = castParams(layer, type)
            %Cast the parameters of the layer.
            % INPUT type: Function handle which is applied to the
            % parameters.
        end

        function params = param2Vec(layer)
            params = [];
        end

        function layer = vec2Param(layer, vec)
        end
    end

end
