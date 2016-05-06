classdef dropoutLayer < Codat.Net.Layer.baseLayer
    %CORRUPTIONLAYER Layer which sets a random fraction of the input array to
    % zero. The dropout mask is generated for the first three dimensions and
    % repeated along the fourth dimension.
    %
    % PROPERTIES
    % p: Fraction of inputs which are set to zero.
    % mask: Mask used for dropout during forward pass. False locations in mask
    %       are dropped.
    % Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

    properties
        p
        mask
    end

    methods
        function obj = dropoutLayer(p)
            obj.p = p;
            obj.isInPlace = true;
        end

        function [ layer, X ] = forward(layer, X)
            layer.mask = rand(size(X),'like',X) < (1 - layer.p);
            X = bsxfun(@times, X, layer.mask)./(1 - layer.p);
        end

        function [DEDX, DEDW] = backward(layer, ~, DEDY)
            DEDX = bsxfun(@times,DEDY,layer.mask)./(1 - layer.p);
            DEDW = [];
        end
    end

end
