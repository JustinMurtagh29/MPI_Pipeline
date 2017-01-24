classdef LocalPol < SynEM.Feature.TextureFeature
    %LOCALPOL Local polynomial filter.
    % PROPERTIES
    % nhood: [Nx1] array of integer specifying the size of the filter in
    %   each dimension.
    % sigma: (Optional) [Nx1] array of float specifying the standard
    %   deviation in each dimension for prior smoothing
    %   (Default: no prior smoothing)
    %
    % Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>
    
    properties
        nhood
        sigma = [];
    end
    
    methods
        function obj = LocalPol(nhood, sigma)
            obj.name = 'LPol';
            obj.nhood = nhood;
            if exist('sigma','var')
                obj.sigma = sigma;
                obj.border = (nhood - 1) + 2.*ceil(3.*sigma);
            else
                obj.border = (nhood - 1);
            end
            obj.numChannels = 1;
            
        end
        
        function feat = calc(obj, raw)
            feat = obj.calculate(raw, obj.nhood, obj.sigma);
        end
    end
    
    methods (Static)
        function feat = calculate(raw, nhood, sigma)
            if ~isempty(sigma)
                raw = SynEM.Feature.GaussFilter.calculate(raw, sigma, ...
                    3, 0, 'same');
            end
            h = ones(nhood)./prod(nhood);
%             raw = raw./255; %legacy, could be removed at some point
            feat = imfilter(raw.*raw,h) - imfilter(raw,h).^2;
        end
    end
    
end

