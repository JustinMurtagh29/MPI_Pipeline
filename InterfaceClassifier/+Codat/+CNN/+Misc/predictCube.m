function [ pred ] = predictCube( cnet, raw, options )
%PREDICTCUBE Predict for input cube.
% INPUT cnet: A Codat.CNN object.
%       raw: 3D input array.
%       options: Options struct with field
%                gpuDev: see cnet.train
%                target_size: 3x1 array.
%                             Input cube is tiled such that each forward
%                             pass produces an output of target_size.
%                             Targets for one cube are stiched together
%                             afterwards. (Default is the whole ROI).
%                convAlg: (Default is current convAlg of cnet)
%                display: (Optional) Positive integer specifying the
%                         printout frequency to track how many tiles have
%                         been processed so far.
%                         (Default: No printout).
%               rotPred: (Optional) Logical flag indicating whether
%                   rotPredict is used for prediction.
%                   (Default: false)
% OUTPUT pred: Prediction for raw.
%
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if options.gpuDev
    cnet = cnet.setParamsToActvtClass(@gpuArray);
end
if ~isfield(options,'target_size') || isempty(options.target_size)
    options.target_size = size(raw);
end
if isfield(options,'convAlg')
    cnet = cnet.setConvMode(options.convAlg);
end
if ~isfield(options,'rotPred')
    options.rotPred = false;
end

%make tiles
[y,x,z] = meshgrid(1:options.target_size(2):size(raw,2) - cnet.border(2), ...
                   1:options.target_size(1):size(raw,1) - cnet.border(1), ...
                   1:options.target_size(3):size(raw,3) - cnet.border(3));
%predict each tile
pred = cell(size(x));
if isfield(options,'display')
    fprintf('[%s] Starting prediction for %d tiles.\n',datestr(now),numel(x));
end
for tile = 1:numel(x)
    currTile = raw(x(tile):min(x(tile) + cnet.border(1) + options.target_size(1) - 1,size(raw,1)), ...
                   y(tile):min(y(tile) + cnet.border(2) + options.target_size(2) - 1,size(raw,2)), ...
                   z(tile):min(z(tile) + cnet.border(3) + options.target_size(3) - 1,size(raw,3)));
    if options.gpuDev
        currTile = gpuArray(currTile);
    end
    if options.rotPred
        pred{tile} = gather(cnet.rotPredict(currTile));
    else
        pred{tile} = gather(cnet.predict(currTile));
    end
    if isfield(options,'display')
        if mod(tile,options.display) == 0
            fprintf('[%s] Finished tile %d out of %d.\n',datestr(now), ...
                     tile,numel(x));
        end
    end
end
pred = cell2mat(pred);

end

