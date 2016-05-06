classdef batchNormalizationLayer < Codat.Net.Layer.baseLayer
    %BATCHNORMALIZATIONLAYER Layer which performs batch normalization.
    % Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

    properties
        sigma2_inf;             %variance for inference
        mu_inf;                 %mean for inference
        gamma;
        beta;
        epsilon = 1e-8;
        muB;                    %mean of current batch
        sigmaB2;                %variance of current batch
        num_fms;
        keepMovingAverage = true;%bool specifying whether an exponential
                                %moving average of mu and sigma is kept
                                %during training
        alpha = 0.5;            %factor for exponential moving average
    end

    methods
        function obj = batchNormalizationLayer(num_fms)
            obj.isInPlace = false;
            obj.numParams = 2*num_fms;
            obj.num_fms = num_fms;
            obj.gamma = ones(1,1,1,num_fms);
            obj.beta = zeros(1,1,1,num_fms);
            obj.mu_inf = zeros(1,1,1,num_fms);
            obj.sigma2_inf = zeros(1,1,1,num_fms);
            obj.numParams = 2*num_fms;
        end

        function [ layer, X ] = forward(layer, X)
            sX = size(X);
            m = prod(sX(1:3));
            layer.muB = layer.sum(X,1:3)./m;
            layer.sigmaB2 = layer.sum(bsxfun(@minus,X,layer.muB).^2,1:3)./m;
            X = bsxfun(@minus,X,layer.muB);
            X = bsxfun(@rdivide,X,sqrt(layer.sigmaB2 + layer.epsilon));
            X = bsxfun(@times,X,layer.gamma);
            X = bsxfun(@plus,X,layer.beta);
            if layer.keepMovingAverage
                layer.mu_inf = layer.alpha.*layer.mu_inf + (1-layer.alpha).*layer.muB;
                layer.sigma2_inf = layer.alpha.*layer.sigma2_inf + (1-layer.alpha).*layer.sigmaB2;
            end
        end

        function [DEDX, DEDW] = backward(layer, X, DEDY)
            sX = size(X);
            m = prod(sX(1:3));
            X = bsxfun(@minus,X,layer.muB);
            DEDXh = bsxfun(@times,DEDY,layer.gamma);
            DEDsigB_2 = -0.5.*layer.sum(DEDXh.*X,1:3).*(layer.sigmaB2 + layer.epsilon).^(-3/2);
            DEDmuB = -layer.sum(DEDXh,1:3)./(sqrt(layer.sigmaB2 + layer.epsilon)) - 2.*DEDsigB_2./m.*layer.sum(X,1:3);
            DEDX = bsxfun(@rdivide,DEDXh,sqrt(layer.sigmaB2 + layer.epsilon)) + bsxfun(@plus,2.*bsxfun(@times,X,DEDsigB_2)./m,DEDmuB./m);
            DEDgamma = layer.sum(DEDY.*bsxfun(@rdivide,X,sqrt(layer.sigmaB2 + layer.epsilon)),1:3);
            DEDbeta = layer.sum(DEDY,1:3);
            DEDW = [DEDbeta(:);DEDgamma(:)];
        end
        
        function vec = param2Vec(layer)
            vec = [layer.beta(:);layer.gamma(:)];
        end
        
        function layer = vec2Param(layer, vec)
            layer.beta = reshape(vec(1:end/2),[1 1 1 layer.num_fms]);
            layer.gamma = reshape(vec((end/2 + 1):end),[1 1 1 layer.num_fms]);
        end
        
        function layer = castParams(layer, type)
            layer.muB = type(gather(layer.muB));
            layer.sigmaB2 = type(gather(layer.sigmaB2));
            layer.mu_inf = type(gather(layer.mu_inf));
            layer.sigma2_inf = type(gather(layer.sigma2_inf));
        end
    end

    methods (Static)
        function X = sum(X,dims)
            for dim = dims
                X = sum(X,dim);
            end
        end
    end

end
