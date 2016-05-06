classdef rmsProp < Codat.Optimizer.optimizer
    %RMSPROP RMSProp optimizer.
    % Rmsprop uses an average of past gradient to normalize the gradient
    % This normalization is a moving average over the root mean squared
    % gradients (T. Tieleman and G. Hinton. Lecture 6.5 - rmsprop: Divide
    % the gradient by a running average of its recent magnitude, 2012.)
    %
    % Properties learningRate: Learning rate by which the gradient is multiplied.
    %            rho: Fraction of how much past information to keep.
    %            momentum: Momentum of optimizer.
    %            epsilon: Conditioning parameter.
    %
    % Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

    properties
        learningRate = 1e-2;
        rho = 0.9;
        epsilon = 1e-6;
        momentum = 0.9;
        eg
        v
    end

    methods
        function rmsOpt = rmsProp(learningRate,rho,momentum,epsilon)
            rmsOpt = rmsOpt@Codat.Optimizer.optimizer;
            rmsOpt.learningRate = learningRate;
            if exist('rho','var') && ~isempty(rho)
                rmsOpt.rho = rho;
            end
            if exist('epsilon','var') && ~isempty(epsilon)
                rmsOpt.epsilon = epsilon;
            end
            if exist('momentum','var') && ~isempty(momentum)
                rmsOpt.momentum = momentum;
            end
        end

        function rmsOpt = init(rmsOpt, numParams)
            rmsOpt.eg = zeros(numParams,1,'single');
            rmsOpt.v = zeros(numParams,1,'single');
        end

        function [paramsNew,rmsOpt] = optimize(rmsOpt,paramsOld,gradient)
            rmsOpt.eg = rmsOpt.rho.*rmsOpt.eg + (1 - rmsOpt.rho)*gradient.^2;
            rmsOpt.v = rmsOpt.momentum.*rmsOpt.v + rmsOpt.learningRate./sqrt(rmsOpt.eg + rmsOpt.epsilon).*gradient;
            paramsNew = paramsOld - rmsOpt.v;
        end

        function optimizer = setParamsToActvtClass(optimizer,actvtClass)
            optimizer.v = actvtClass(gather(optimizer.v));
            optimizer.eg = actvtClass(gather(optimizer.eg));
        end
    end

end
