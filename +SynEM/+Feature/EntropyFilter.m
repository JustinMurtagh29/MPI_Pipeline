classdef EntropyFilter < SynEM.Feature.TextureFeature
    %ENTROPYFILTER Local entropy filter.
    % PROPERTIES
    % nhood: [Nx1] array of integer specifying the size of the filter in
    %   each dimension.
    % sigma: (Optional) [Nx1] array of float specifying the standard
    %   deviation in each dimension for prior smoothing
    %   (Default: no prior smoothing)
    % normalize: (Optional) Logical indicating whether raw should be
    %   normalized to [0, 255].
    %   (Default: false)
    %
    % NOTE raw is cast to uint8 before calculating entropyfilt which might make
    %      the normalization important.
    %
    % Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>
    
    properties
        nhood
        sigma = [];
        normalize = false;
    end
    
    methods
        function obj = EntropyFilter(nhood, sigma, normalize)
            obj.name = 'Entropy';
            obj.nhood = nhood;
            if exist('sigma','var') && ~isempty(sigma)
                obj.sigma = sigma;
                obj.border = (nhood - 1) + 2.*ceil(3.*sigma);
            else
                obj.border = (nhood - 1);
            end
            if exist('normalize','var') && ~isempty(normalize)
                obj.normalize = normalize;
            end
            obj.numChannels = 1;
            
        end
        
        function feat = calc(obj, raw)
            feat = obj.calculate(raw, obj.nhood, obj.sigma, obj.normalize);
        end
    end
    
    methods (Static)
        function feat = calculate(raw, nhood, sigma, normalize)
            if ~isempty(sigma)
                raw = SynEM.Feature.GaussFilter.calculate(raw, sigma, 3,...
                    0, 'same');
            end
            if normalize
                raw = ((raw - min(raw(:)))./ ...
                    (max(raw(:)) - min(raw(:)))).*255;
            end
            feat = cast(entropyfilt(uint8(raw),ones(nhood)),'like',raw);
        end
    end
    
end

