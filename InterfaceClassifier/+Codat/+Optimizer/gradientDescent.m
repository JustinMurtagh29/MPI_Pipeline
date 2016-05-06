classdef gradientDescent < Codat.Optimizer.optimizer
    %GRADIENTDESCENT Optimization by gradient descent.
    % Properties learningRate: Prefactor of gradient.
    %            momentum: Momentum strength
    % Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>
    
    properties
        learningRate
        momentum
        v
    end
    
    methods
        function gd = gradientDescent(learningRate,momentum)
            gd = gd@Codat.Optimizer.optimizer;
            gd.learningRate = learningRate;
            gd.momentum = momentum;
        end
        
        function gd = init(gd, numParams)
            gd.v = zeros(numParams,1,'single');
        end
        
        function [paramsNew, gd] = optimize(gd,paramsOld,gradient)
            gd.v = gd.learningRate.*gradient + gd.momentum.*gd.v;
            paramsNew = paramsOld - gd.v;
        end
        
        function optimizer = setParamsToActvtClass(optimizer,actvtClass)
            optimizer.v = actvtClass(gather(optimizer.v));
        end
    end
    
end

