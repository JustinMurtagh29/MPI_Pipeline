function [ cnet, W, b ] = convertSegEM( segEM )
%CONVERTSEGEM Convert SegEM CNN to Codat CNN.
% INPUT segEM: SegEM cnn object.
% OUTPUT cnet: Codat cnn object with parameters from segEM.
%        W,b: Parameters of segEM in format for Codat cnn.
% NOTE Check optimizer settings before continuing to train the Codat.CNN.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

cLayer = segEM.numLayer - 1;
featureMaps = [1 segEM.numFeature 1];
filterSize = {segEM.filterSize};
dropout = 0;
stride = {[1 1 1]};
maxPool = false;
shortcut = zeros(1,segEM.numLayer);
batchNorm = false(1,segEM.numLayer);
optimizer = Codat.Optimizer.gradientDescent(1e-8,0.9);
cnet = Codat.CNN.cnn(cLayer, featureMaps, filterSize, stride, maxPool, dropout, shortcut, batchNorm, 'tanh', 'squared', optimizer);
W = cell(1,cnet.layer);
b = cell(1,cnet.layer);
for lyr = 2:cnet.layer
    W{lyr} = zeros(size(cnet.W{lyr}),'like',cnet.W{lyr});
    for prevFm = 1:cnet.featureMaps(lyr - 1)
        for fm = 1:cnet.featureMaps(lyr)
            W{lyr}(:,:,:,prevFm,fm) = segEM.layer{lyr}.W{prevFm,fm};
        end
        b{lyr} = segEM.layer{lyr}.B';
    end
end
W = cellfun(@(x)flip(x,4),W,'UniformOutput',false);
cnet.W = W;
cnet.b = b;




end

