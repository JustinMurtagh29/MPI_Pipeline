function [ ac, ac_ref ] = synSpatialFreq( raw, raw_ref, structureSize, ...
    voxelSize, mode )
%SYNSPATIALFREQ
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if ~exist('voxelSize', 'var') || isempty(voxelSize)
    voxelSize = [11.24, 11.24, 28];
end
if ~exist('mode', 'var') || isempty(mode)
    mode = 'gaussbpf';
end

% fft transform and shift spectrum
fraw = (fftn(single(raw)));
fraw_ref = (fftn(single(raw_ref)));

% autocorr
acRaw = fftshift(ifftn(fraw.*conj(fraw)));
acRaw_ref = fftshift(ifftn(fraw_ref.*conj(fraw_ref)));

% normalize
acRaw = acRaw./max(acRaw(:));
acRaw_ref = acRaw_ref./max(acRaw_ref(:));

% filter raw
switch mode
    case 'gaussbpf'
        filter = gaussBpf(size(acRaw), structureSize, voxelSize);
    case 'hardbpf'
        filter = hardBpf(size(acRaw), structureSize, voxelSize);
    case 'pixelbpf'
        filter = pixelBpf(size(acRaw), structureSize, voxelSize);
end
acRaw = acRaw.*filter;

% filter raw_ref
switch mode
    case 'gaussbpf'
        filter = gaussBpf(size(acRaw_ref), structureSize, voxelSize);
    case 'hardbpf'
        filter = hardBpf(size(acRaw_ref), structureSize, voxelSize);
    case 'pixelbpf'
        filter = pixelBpf(size(acRaw_ref), structureSize, voxelSize);
end
acRaw_ref = acRaw_ref.*filter;

% mean
ac = mean(acRaw(acRaw > 0));
ac_ref = mean(acRaw_ref(acRaw_ref > 0));
end


function filter = pixelBpf(siz, structureSize, vxS)
[x,y,z] = ndgrid(0:(siz(1) - 1), 0:(siz(2) - 1), 0:(siz(3) - 1));
x = (x - floor(siz(1)/2)).*vxS(1);
y = (y - floor(siz(2)/2)).*vxS(2);
z = (z - floor(siz(3)/2)).*vxS(3);
h1 = x.^2 + y.^2 + z.^2 >= structureSize(1)^2;
h2 = x.^2 + y.^2 + z.^2 <= structureSize(2)^2;
filter = single(h1.*h2);

end


function [sigL, sigH] = determineStd(siz, structureSize, voxelSize)
% wavelength of structures (structure size is wavelength/2)
lowerWV = 2.*bsxfun(@rdivide,structureSize(:,1),voxelSize);
upperWV = 2.*bsxfun(@rdivide,structureSize(:,2),voxelSize);

%pixel frequency of structures
lowerF = 1./upperWV;
upperF = 1./lowerWV;

%frequency increment per pixel in fourier space
dF = 1./siz;

%region in fft space
sigH = bsxfun(@rdivide,lowerF,dF);
sigL = bsxfun(@rdivide,upperF,dF);

end

function filter = hardBpf(siz, structureSize, voxelSize)

% frequency increment per pixel in fourier space
dF = 1./voxelSize./siz; % sampling frequency / length (in 1/nm)

% frequency of structure
maxF = 1./(2.*structureSize(1)); % in 1/nm
minF = 1./(2.*structureSize(2)); % in 1/nm

%define grid for filter
[x,y,z] = ndgrid(0:(siz(1) - 1), 0:(siz(2) - 1), 0:(siz(3) - 1));
x = abs(x - floor(siz(1)/2)).*dF(1);
y = abs(y - floor(siz(2)/2)).*dF(2);
z = abs(z - floor(siz(3)/2)).*dF(3);

%define filter
highPass = single(x.^2 + y.^2 + z.^2 >= minF.^2);
lowPass = single(x.^2 + y.^2 + z.^2 <= maxF.^2);
filter = lowPass.*highPass;
end

function filter = gaussBpf(siz, structureSize, voxelSize)

[sigL, sigH] = determineStd(siz, structureSize, voxelSize);

%define grid for filter
[x,y,z] = ndgrid(0:(siz(1) - 1), 0:(siz(2) - 1), 0:(siz(3) - 1));
x = (x - floor(siz(1)/2)).^2;
y = (y - floor(siz(2)/2)).^2;
z = (z - floor(siz(3)/2)).^2;

%prepare sig
sigH = sigH.^2;
sigL = 2.*sigL.^2;

%define filter
lowPass = exp(-x./sigL(1) - y./sigL(2) - z./sigL(3));
highPass = 1 - exp(-x./sigH(1) - y./sigH(2) - z./sigH(3));
filter = lowPass.*highPass;
end