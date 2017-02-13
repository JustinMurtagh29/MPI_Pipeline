classdef AverageFilter < SynEM.Feature.TextureFeature
    %AVERAGEFILTER Local average filter
    %
    % PROPERTIES
    % type: String
    %   Type specification. Options are
    %       'cube'
    %       'ball'
    % parameters: Cell array containing type specific paramters.
    %   'cube' - edgeLength: [Nx1] array of integer specifying the cube
    %       length in each dimension.
    %   'ball' - radius: Float specifying the ball radius.
    %          - voxelSize: (Optional) [Nx1] array of float specifying the
    %               voxel size in each dimension.
    %               (Default: 1)
    % Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>
    
    properties
        type
        parameters
    end
    
    methods
        function obj = AverageFilter(type, varargin)
            obj.name = 'AvgFilt';
            obj.type = type;
            switch type
                case 'cube'
                    obj.parameters{1} = varargin{1};
                    obj.border = (varargin{1} - 1);
                case 'ball'
                    obj.parameters = varargin;
                    r = varargin{1};
                    if length(varargin) == 2
                        sc = varargin{2};
                    else
                        sc = ones(3,1);
                    end
                    fSize = ceil(r./sc);
                    obj.border = 2*fSize + 1;
                otherwise
                    error('Unknown type %s.', type);
            end
            obj.numChannels = 1;
            if iscolumn(obj.border)
                obj.border = obj.border';
            end
        end
        
        function feat = calc(obj, raw)
            feat = obj.calculate(raw, obj.type, obj.parameters{:});
        end
    end
    
    methods (Static)
        function feat = calculate(raw, type, varargin)
            switch type
                case 'cube'
                    nD = ndims(raw);
                    l = varargin{1};
                    if length(l) == 1
                        l = repmat(l,ones(nD,1));
                    end
                    pL = prod(l);

                    feat = raw;
                    for dim = 1:nD
                        sz = ones(1, nD);
                        sz(dim) = l(dim);
                        h = reshape(ones(l(dim),1),sz);
                        feat = convn(feat, h, 'same');
                    end
                    feat = feat./pL;
                case 'ball'
                    r = varargin{1};
                    if length(varargin) == 2
                        sc = varargin{2};
                    else
                        sc = ones(3,1);
                    end
                    fSize = ceil(r./sc);
                    [x,y,z] = meshgrid(-fSize(1):fSize(1), ...
                        -fSize(2):fSize(2), -fSize(3):fSize(3));
                    h = (sc(1).*x).^2 + (sc(2).*y).^2 + (sc(3).*z).^2 ...
                        <= r^2;
                    h = h/sum(h(:));
                    feat = imfilter(raw,double(h));
                otherwise
                    error('Type %s is not defined.', type);
            end
        end
    end
    
end

