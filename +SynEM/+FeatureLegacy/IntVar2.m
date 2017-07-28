classdef IntVar2 < SynEM.Feature.TextureFeature
    % IntVar Intensity/variance filter.
    % PROPERTIES
    % nhood: [Nx1] array of integer specifying the size of the filter in
    %   each dimension.
    % sigma: (Optional) [Nx1] array of float specifying the standard
    %   deviation in each dimension for prior smoothing
    %   (Default: no prior smoothing)
    % n_mean: The mean used for raw data normalization.
    % n_std: The standard deviation used for raw data normlaization.
    %
    % NOTE n_mean and n_std are used to ensure a consistent normalization
    %      w.r.t. to uint8 raw data.
    %
    % Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>
    
    properties
        nhood
        n_mean = 0
        n_std = 1
    end
    
    methods
        function obj = IntVar2(nhood, n_mean, n_std)
            obj.name = 'IntVar';
            obj.nhood = nhood;
            obj.border = (nhood - 1);
            obj.numChannels = 1;
            if exist('n_mean', 'var') && ~isempty(n_mean)
                obj.n_mean = n_mean;
            end
            if exist('n_std', 'var') && ~isempty(n_std)
                obj.n_std = n_std;
            end
        end
        
        function feat = calc(obj, raw)
            feat = obj.calculate(raw, obj.nhood, obj.n_mean, obj.n_std);
        end
    end
    
    methods (Static)
        function feat = calculate(raw, nhood, n_mean, n_std)
            h = ones(nhood, 'like', raw);
            raw = (raw.*n_std + n_mean)./255;
            feat = convn(raw.*raw, h, 'same') - convn(raw, h, 'same').^2;
        end
    end
    
end

