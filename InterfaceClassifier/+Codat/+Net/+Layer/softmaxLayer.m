classdef softmaxLayer < Codat.Net.Layer.baseLayer
    %SIGMOIDLAYER Softmax nonlinearity along fourth dimension of input.
    %
    % PROPERTIES
    % onlyBackpass: Bool specifying whether error should be propagated back
    %       through non-linarity or simply passed back to the next layer. This
    %       options can be used if the softmax layer is followed by a softmax
    %       loss layer to simplify computations.
    %
    % NOTE This layer is only intended to be used directly before a sofmax
    % loss layer with onlyBackpass set to true. The general backpropagation
    % is implemented yet very inefficient.
    %
    % Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

    properties
        onlyBackpass = true
    end

    methods
        function obj = softmaxLayer(onlyBackpass)
            if exist('onlyBackpass','var') || ~isempty(onlyBackpass)
                obj.onlyBackpass = onlyBackpass;
            end
            obj.isInPlace = true;
        end

        function [ layer, Y ] = forward(layer, X)
            num = exp(X);
            denom = sum(num,4);
            Y = bsxfun(@rdivide,num,denom);
        end

        function [DEDX, DEDW] = backward(layer, Y, DEDY)
            if layer.onlyBackpass
                DEDX = DEDY;
            else
                %TODO: optimize if ever needed
                DEDX = zeros(size(DEDY),'like',DEDY);
                for x = 1:size(Y,1)
                    for y = 1:size(Y,2)
                        for z = 1:size(Y,3)
                            J = - permute(Y(x,y,z,:),[4 1 2 3])*permute(Y(x,y,z,:),[1 4 2 3]) + diag(permute(Y(x,y,z,:),[4 1 2 3]));
                            DEDX(x,y,z,:) = permute(J*permute(DEDY(x,y,z,:),[4 1 2 3]),[4 3 2 1]);
                        end
                    end
                end
            end
            DEDW = [];
        end
    end

end
