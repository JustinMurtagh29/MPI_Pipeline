function [ raw, target, targetWeights ] = loadDataRMMV( path, weight, sT )
%LOADDATARMMV Load for membrane prediction with RMMV input (raw, membrane,
%mitos, vesicles).
% INPUT path: Path to training data file.
%       weight: (Optional) Double weight factor for membran class.
%               (Default: Emprical weight to balance classes)
%       sT: (Optional) Logical specifying whether to train on small target.
%           (Default: false)
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

m = load(path);

if ~exist('weight','var') || isempty(weight)
    weight = sum(m.target(:) == 1)/sum(m.target(:) == -1);
end

%load road and mmv
raw = (single(m.raw) - 122)./22;
mmv = m.feat;

%crop raw to mmv and cat both together
sR = size(raw);
sF = size(mmv);
border = (sR - sF(1:3));
raw = raw(border(1)/2 + 1:end - border(1)/2, ...
              border(2)/2 + 1:end - border(2)/2, ...
              border(3)/2 + 1:end - border(3)/2);
raw = cat(4,raw,mmv);

target = single(m.target);
if exist('sT','var') && sT
    target = imdilate(target,ones(2,2));
end
targetWeights = ones(size(target),'single');
targetWeights(target == -1) = weight;
targetWeights(target == 0) = 0;


end

