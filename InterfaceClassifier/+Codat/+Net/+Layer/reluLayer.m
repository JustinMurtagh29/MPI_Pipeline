classdef reluLayer < Codat.Net.Layer.baseLayer
    %RELULAYER Pointwise rectified linear unit.
    % Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>
    
    properties
    end
    
    methods
        function obj = reluLayer()
            obj.isInPlace = true;
        end
        
        function [layer, X] = forward(layer, X)
            X = cast(max(0,X),'like',X); %is cast necessary here?
        end
        
        function [DEDX, DEDW] = backward(~, Y, DEDY)
            DEDX = DEDY.*cast(Y > 0,'like',Y);
            DEDW = [];
        end
    end
    
end

