function preds = predictStacks( cnet, stacks, options )
%PREDICTSTACKS Predict multiple stacks
% INPUT cnet: A Codat.CNN object.
%       stacks: Cell array with paths to stacks.
%       options: Options struct with field
%                gpuDev: see cnet.train
%                data_pre: see cnet.train
%                target_size: 3x1 array.
%                             Input cube is tiled such that each forward
%                             pass produces an output of target_size.
%                             Targets for one cube are stiched together
%                             afterwards. (Default is the whole stack).
%                val_fwd_alg (Default is current fwdAlg of cnet)
% OUTPUT preds: Cell array of predictions for corresponding stack.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

%set options
if ischar(options.data_pre) %in case validate is called directly
    options.data_pre = str2func(options.data_pre);
end
if options.gpuDev
%    gpuDevice(options.gpuDev);
    cnet = cnet.setParamsToActvtClass(@gpuArray);
end
if ~isfield(options,'target_size') || isempty(options.target_size)
    autoSize = true;
else
    autoSize = false;
end
if isfield(options,'val_fwd_alg')
    cnet = cnet.setConvMode(options.val_fwd_alg);
end

%predict stacks
preds = cell(length(stacks),1);
for i = 1:length(stacks)
    fprintf('[%s] Prediction for stack %d/%d.\n',datestr(now),i,length(stacks));
    raw = options.data_pre(stacks{i});
    if autoSize
        options.target_size = size(raw) - cnet.border;
    end
    if options.gpuDev
        raw = gpuArray(raw);
    end
    
    tilesPerDim = ceil((size(raw) - cnet.border)./options.target_size);
    predTiles = cell(tilesPerDim);
    numTiles = prod(tilesPerDim);
    for tile = 1:numTiles
        cubeTile = tileCube(raw,options.target_size,cnet.border,tile);
        predTiles{tile} = gather(cnet.predict(cubeTile));
    end
    preds{i} = cell2mat(predTiles);
end


end

function cubeTile = tileCube(cube,targetSize,border,tileNo)
cubeSize = size(cube);
%define (1,1,1) coordinate of each cube tile
[y,x,z] = meshgrid(1:targetSize(1):cubeSize(1) - border(1), ...
                   1:targetSize(2):cubeSize(2) - border(2), ...
                   1:targetSize(3):cubeSize(3) - border(3));
cubeTile = cube(x(tileNo):min((x(tileNo) + border(1) + targetSize(1) - 1),cubeSize(1)), ...
                y(tileNo):min((y(tileNo) + border(2) + targetSize(2) - 1),cubeSize(2)), ...
                z(tileNo):min((z(tileNo) + border(3) + targetSize(3) - 1),cubeSize(3)),:);

end

