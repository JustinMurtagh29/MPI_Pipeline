function [ coords, az, el ] = randomDirection( seedPoint, l, sampleDist )
%RANDOMDIRECTION Draw a line in a random direction from a seed point.
% INPUT seedPoint: [1x3] int
%           Cartesian coordinates of the seed point.
%       l: double
%           Length of the line.
%       sampleDist: (Optional) double
%           Distance of the sampled points on the line.
%           (Default: l/10)
% OUTPUT coords: [Nx3] int
%           Points sampled from the line.
%        az: double
%           The azimuth angle of the line.
%        el: double
%           The elevation angle of the line.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

az = rand(1)*2*pi;
el = rand(1)*pi - pi/2;

if ~exist('sampleDist', 'var') || isempty(sampleDist)
    sampleDist = l/10;
end

[x, y, z] = sph2cart(az, el, (0:sampleDist:l)');
coords = [x, y, z];
coords = round(bsxfun(@plus, coords, seedPoint(:)'));

end

