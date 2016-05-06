function [x_train,y_train, mask] = tileTrainingCubes( inputCube,labelCube,targetSize,border,maxTiles,tileNumbers, targetMask )
%TILETRAININGCUBES Divide labelCube into tiles of targetSize and get the
%corresponding raw input.
% INPUT inputCube: 4D input cube. (4th dimension corresponds to input
%                  channels).
%       targetCube: 4D cube with target labels centered in inputCube (4th
%                   dimension corresponds to output channels).
%       targetSize: Size of the desired targets for CNNMPF
%       border: The CNNMPF border property
%       maxTiles: Integer specifying maximal number of tiles to produce.
%       tileNumber: Array specifying tiles to extract. tileNumber can range
%           from 1 to max(totalNumberTiles,maxTiles).
%       targetMask: Mask for target values.
% OUTPUT: x_train: Cell-array of input training cubes.
%         y_train: Cell-array of target training cubes.
%         mask: Cell-array of targetMask cubes.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

%calculate how often a targetSize cube fits into the targetCube
targetCubeSize = size(labelCube);
noTiles = floor(targetCubeSize(1:3)./targetSize);
noTrainingPoints = prod(noTiles);
if ~exist('maxTiles','var') || isempty(maxTiles)
    maxTiles = noTrainingPoints;
else
    maxTiles = min(maxTiles,noTrainingPoints);
end
if ~exist('tileNumbers','var') || isempty(tileNumbers)
    tileNumbers = 1:maxTiles;
elseif ~isrow(tileNumbers)
    tileNumbers = tileNumbers';
end
if ~exist('targetMask','var') || isempty(targetMask)
    if nargout == 3
        error('Too many output arguments. Specify a targetMask for mask output.')
    end
end

x_train = cell(length(tileNumbers),1);
y_train = cell(length(tileNumbers),1);
mask = cell(length(tileNumbers),1);
coordT = (size(inputCube(:,:,:,1)) - targetCubeSize(1:3))./2;
border = border./2;
if any(border > coordT)
    error('Specified border of [%s] is too large. The maximum possible border for this training data is [%s].',num2str(2.*border),num2str(2.*coordT));
end

%cut out the small cubes and the corresponding input cubes
[x,y,z] = meshgrid(1:targetSize(1):noTiles(1)*targetSize(1),1:targetSize(2):noTiles(2)*targetSize(2),1:targetSize(3):noTiles(3)*targetSize(3));

for i = 1:length(tileNumbers)
    tile = tileNumbers(i);
    y_train{i} = labelCube(x(tile) : (x(tile) + targetSize(1) - 1), ...
                           y(tile) : (y(tile) + targetSize(2) - 1), ...
                           z(tile) : (z(tile) + targetSize(3) - 1),:);
    if exist('targetMask','var') && ~isempty(targetMask)
        mask{i} = targetMask(x(tile) : (x(tile) + targetSize(1) - 1), ...
                             y(tile) : (y(tile) + targetSize(2) - 1), ...
                             z(tile) : (z(tile) + targetSize(3) - 1),:);
    end
    x_train{i} = inputCube(x(tile) + coordT(1) - border(1) : x(tile) + coordT(1) + targetSize(1) + border(1) - 1, ...
                           y(tile) + coordT(2) - border(2) : y(tile) + coordT(2) + targetSize(2) + border(2) - 1, ...
                           z(tile) + coordT(3) - border(3) : z(tile) + coordT(3) + targetSize(3) + border(3) - 1,:);
    if tile > maxTiles
        return;
    end
end

end

