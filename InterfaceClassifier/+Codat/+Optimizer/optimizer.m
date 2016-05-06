classdef optimizer
    %OPTIMIZER Base class for an optimizer.
    % Each optimizer is intended to receive a vector of all parameters in
    % the cnn class and a gradient vector and output a vector of updated
    % parameters.
    % Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>
    
    properties
    end
    
    methods (Abstract)
        [paramsNew,optimizer] = optimize(optimizer,paramsOld,gradient);
        optimizer = setParamsToActvtClass(optimizer,actvtClass);
        optimizer = init(optimizer, numParams);
    end
    
end

