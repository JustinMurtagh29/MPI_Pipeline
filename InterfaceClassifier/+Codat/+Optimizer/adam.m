classdef adam < Codat.Optimizer.optimizer
    %ADAM Implementation of Adam (Kingma, D., & Ba, J. (2014). Adam: A
    %method for stochastic optimization. arXiv preprint arXiv:1412.6980.)
    % Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>
    
    properties
        t = 0;
        stepsize = 0.0002;
        decay_mom1 = 0.9;
        decay_mom2 = 0.999;
        offset = 1e-8;
        lambda = 1 - 1e-8;
        est_mom1_b
        est_mom2_b
        learningRate = 1
    end
    
    methods
        function ad = adam(stepsize,lambda,decay_mom1,decay_mom2,offset)
            ad = ad@Codat.Optimizer.optimizer;
            if exist('stepsize','var') && ~isempty(stepsize)
                ad.stepsize = stepsize;
            end
            if exist('lambda','var') && ~isempty(lambda)
                ad.lambda = lambda;
            end
            if exist('decay_mom1','var') && ~isempty(decay_mom1)
                ad.decay_mom1 = decay_mom1;
            end
            if exist('decay_mom2','var') && ~isempty(decay_mom2)
                ad.decay_mom2 = decay_mom2;
            end
            if exist('offset','var') && ~isempty(offset)
                ad.offset = offset;
            end
        end
        
        function ad = init(ad, numParams)
            ad.est_mom1_b = zeros(numParams,1,'single');
            ad.est_mom2_b = zeros(numParams,1,'single');
        end
        
        function [paramsNew,ad] = optimize(ad,paramsOld,gradient)
            ad.t = ad.t + 1;
            decay_mom1_t = ad.decay_mom1*ad.lambda^(ad.t-1);
            ad.est_mom1_b = decay_mom1_t.*ad.est_mom1_b + (1 - decay_mom1_t).*gradient;
            ad.est_mom2_b = ad.decay_mom2.*ad.est_mom2_b + (1 - ad.decay_mom2).*gradient.^2;
            alpha_t = ad.stepsize*sqrt(1 - ad.decay_mom2^ad.t)/(1 - ad.decay_mom1^ad.t);
            paramsNew = paramsOld - ad.learningRate.*alpha_t*ad.est_mom1_b./(sqrt(ad.est_mom2_b + ad.offset));
        end
        
        function optimizer = setParamsToActvtClass(optimizer,actvtClass)
            optimizer.est_mom1_b = actvtClass(gather(optimizer.est_mom1_b));
            optimizer.est_mom2_b = actvtClass(gather(optimizer.est_mom2_b));
        end
    end
    
end

