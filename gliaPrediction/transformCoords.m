function [ coords ] = transformCoords( point,bbox,local )
%TRANSFORMCOORDS Summary of this function goes here
%   Detailed explanation goes here

    if local
        coords = point + bbox(:,1)' - [1 1 1];
    end
    if ~local
        coords = point - bbox(:,1)' + [1 1 1];
    end
end

