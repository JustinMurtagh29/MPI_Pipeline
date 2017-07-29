function [ stackD ] = downsampleStack( stack, raw, sig )
%DOWNSAMPLESTACK Downsample a stack knott data.
% INPUT stack: 3d uint8
%           Knott data (5x5x5 nm resolution). Will be downsampled to
%           11.24 x 11.24 x 28 nm.
%       raw: 3d uint8
%           Sample of our data used for histogram equalization.
%       sig: (Optional) double
%           Standard deviation for gaussian smoothing filter.
%           (Default: no prior smoothing).
% OUTPUT stackD: 3d uint8
%           Downsampled stack.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

stack = single(stack);
if exist('sig', 'var') && ~isempty(sig)
    stack = SynEM.Feature.GaussFilter.calculate(stack, ...
        sig, 3, 0, 'valid');
end

% %downsample to our resolution
% Xq = 1:11.24/5:size(stack, 2);
% Yq = 1:11.24/5:size(stack, 1);
% Zq = 1:28/5:size(stack, 3);
% stackD = zeros(length(Yq), length(Xq), length(Zq));

%to conserve memory do it iteratively
%s = 50;
%for i = 1:ceil(size(Zq, 2)/s)
%    idx = ((i-1)*s + 1):min(i*s,size(Zq, 2));
%    zIdx = max(1, Zq(idx(1)) - 2):min(size(stack, 3), Zq(idx(end)) + 2);
%    stackD(:,:,idx) = interp3(single(stack(:,:,zIdx)), Xq(:), ...
%        Yq(:), Zq(idx) - zIdx(1) + 1, 'cubic');
%end

%direct interpolation (not iteratively)
% stackD = Visualization.interpVid(stack, [5, 5, 5]./[20, 20, 35], ...
%     'linear');
% stackD = Visualization.interpVid(stackD, ...
%     [20, 20, 35]./[11.24, 11.24, 28], 'linear');
stackD = Visualization.interpVid(stack, [5, 5, 5]./[11.24, 11.24, 28], ...
    'linear');

%make psd less pronounced
% minV = 90;
% stackD(stackD < minV) = minV;
softT = 110;
tmp = 1.3.*(softT - stackD(stackD < softT));
stackD(stackD < softT) = stackD(stackD < softT) + ...
    (rand(length(tmp), 1)./2 + 0.5).*tmp;

%some noise
stackD = uint8(stackD + randn(size(stackD),'single')*15); % a bit of noise

if exist('raw', 'var') && ~isempty(raw)
    hgram = histcounts(raw, -0.5:1:255.5);
    stackD = reshape(histeq(stackD(:), hgram), size(stackD));
end

end

