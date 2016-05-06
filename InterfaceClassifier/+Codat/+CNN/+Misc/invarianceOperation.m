function [ raw, target, varargout ] = invarianceOperation( i, raw, target, varargin )
%INVARIANCEOPERATION Apply i-the invariance operation to raw, target and
%additional cubes in varargin. The first four operations rotate each image
%in the x-y-plane by (i-1)*90 degrees
% INPUT i: i-th invariance operation. The first four operations rotate each
%          image in the x-y-plane by (i-1)*90 degrees. The 5-8th operation
%          do the same for the cube flipped along the z-dimension.
%       raw: 3D or 4D input cube.
%       target: 3D or 4D target cube.
%       varargin: Arbitrary number of additional cubes.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

varargout = cell(length(varargin),1);
if i <= 4
    raw = rot90(raw,i-1);
    target = rot90(target,i-1);
    for l = 1:length(varargin)
        varargout{l} = rot90(varargin{l},i-1);
    end
elseif i > 4
    raw = flip(raw,3);
    target = flip(target,3);
    raw = rot90(raw,i-5);
    target = rot90(target,i - 5);
    for l = 1:length(varargin)
        varargin{i} = flip(varargin{l},3);
        varargout{l} = rot90(varargin{l},i - 5);
    end
end


end

