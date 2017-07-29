classdef StdFilter < SynEM.Feature.TextureFeature
    %STDFILTER Local standard deviation filter.
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
        function obj = StdFilter(nhood, sigma)
            obj.name = 'StdFilt';
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
            
            % matlab stdandard filter (internally converts to double and
            % thus has higher ram usage)
%             feat = stdfilt(raw,ones(nhood));
            
            % do standard filtering by hand
            n = prod(nhood);
            conv1 = convn(raw.*raw, ones(nhood(1), 1, 1, 'like', raw), 'same');
            conv1 = conv1./(n - 1);
            conv1 = convn(conv1, ones(1, nhood(2), 1, 'like', raw), 'same');
            conv1 = convn(conv1, ones(1, 1, nhood(3), 'like', raw), 'same');
            conv2 = convn(raw, ones(nhood(1), 1, 1, 'like', raw), 'same');
            conv2 = conv2./sqrt(n*(n - 1));
            conv2 = convn(conv2, ones(1, nhood(2), 1, 'like', raw), 'same');
            conv2 = convn(conv2, ones(1, 1, nhood(3), 'like', raw), 'same');
            feat = sqrt(max((conv1 - conv2.*conv2),0));
        end
    end
end

