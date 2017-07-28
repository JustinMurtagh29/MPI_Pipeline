classdef StructureTensor < SynEM.Feature.TextureFeature
    %STRUCTURETENSOR Multi-dimensional structure tensor.
    %
    % PROPERTIES
    % sigmaD: Integer scalar or vector specifying the standard deviation
    %       of the gaussian derivative for all or seperately for each
    %       dimension of the input image.
    % sigmaW: Integer scalar or vector specifying the standard deviation
    %       of the gaussian window for all or seperately for each
    %       dimension of the input image.
    % truncateD: (Optional) Double. Truncate filter at this many standard
    %         deviations.
    %         (Default: 3)
    % truncateW: (Optional) Double. Truncate filter at this many standard
    %         deviations.
    %         (Default: 3)
    % convMode: String specifying the convolution mode.
    %       (see convn).
    %
    % Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>
    
    properties
        sigmaD
        sigmaW
        truncateD = 3;
        truncateW = 3;
        convMode = 'same'
    end
    
    methods
        function obj = StructureTensor(sigmaD, sigmaW, truncateD, ...
                truncateW, convMode)
            obj.name = 'StructureTensor';
            obj.sigmaD = sigmaD;
            obj.sigmaW = sigmaW;
            if exist('truncateD','var') && ~isempty(truncateD)
                obj.truncateD = truncateD;
            end
            if exist('truncateW','var') && ~isempty(truncateW)
                obj.truncateW = truncateW;
            end
            if exist('convMode','var') && ~isempty(convMode)
                obj.convMode = convMode;
            end
            obj.numChannels = length(sigmaD)*(length(sigmaD) - 1)/2;
            obj.border = 2.*(ceil(sigmaW*obj.truncateW) + ...
                ceil(sigmaD*obj.truncateD));
        end
        
        function S = calc(obj, raw)
            S = obj.calculate(raw, obj.sigmaD, obj.sigmaW, ...
                obj.truncateD, obj.truncateW, obj.convMode);
        end
    end
    
    methods (Static)
        function S = calculate(raw, sigmaD, sigmaW, truncateD, ...
                truncateW, convMode)
            nD = ndims(raw);
            grad = cell(nD,1);
            for dim = 1:ndims(raw)
                order = zeros(nD,1);
                order(dim) = 1;
                grad{dim} = SynEM.Feature.GaussFilter.calculate(raw, ...
                    sigmaD, truncateD, order, convMode);
            end

            S = cell(nD,nD);
            for i = 1:nD
                for j = i:nD
                    S{i,j} = SynEM.Feature.GaussFilter.calculate( ...
                        grad{i}.*grad{j}, sigmaW, truncateW, 0, convMode);
                end
            end
        end
    end
    
end

