function [feat, border, numFeatures] = localAverage( raw, type, varargin )
%LOCALAVERAGE Calculate a multidimensional local average filter.
% INPUT raw: n or 3 -dimensional input array.
%       type: String specifying the shape of the average filter. Options are
%           'cube': Calculate the average over a cube with edge length
%               defined by an integer vector for each dimension or by a
%               scalar for all dimensions in varargin. Cube length should
%               be odd in all dimensions but this is not enforced.
%           'ball': Calculate the average over a sphere with a specified
%               radius as the first varargin argument. You can define the
%               scaling of pixels in each dimension via an optional vector
%               as second input argument (e.g. specify [1, 1, 2.5] if the
%               size of third dimension should be scaled by a factor of 2.5).
%               This filter is only implemented for 3-d input arrays.
% OUTPUT feat: Array of the same size as raw containing the filter response.
%        border: Integer array specifying the total border of the filter
%                in each dimension.
%        numFeatures: Resulting feature space size per voxel.
%
% NOTE It is possible to calculate only the last two outputs by calling the
%      function without the raw input argument.
%
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

switch type
    case 'cube'
        border = (varargin{1} - 1);
        numFeatures = 1;
        if isempty(raw)
            feat = [];
            return;
        end
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
        [x,y,z] = meshgrid(-fSize(1):fSize(1),-fSize(2):fSize(2), ...
            -fSize(3):fSize(3));
        h = (sc(1).*x).^2 + (sc(2).*y).^2 + (sc(3).*z).^2 <= r^2;
        h = h/sum(h(:));
        border = (size(h) - 1);
        numFeatures = 1;
        if isempty(raw)
            feat = [];
            return;
        end
        feat = imfilter(raw,double(h));
    otherwise
        error('Type %s is not defined.', type);
end


end
