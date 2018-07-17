function job = bigFwdPassMyelin( p, bbox )
% job = bigFwdPass( p, bbox )
%   Applies the CNN for myelin detection to the raw data
%   within the specified bounding box.

% This is not of any importance due to CNN translation invariance,
% can be choosen for computational efficency, currenlty optimized for
% running on GPU with 12GB, should be multiples of 128, this is same
% as tileSize right now, no reason it has to be.
% Author: Manuel Berning <manuel.berning@brain.mpg.de>

cubeSize = [512 512 256];
assert(all(mod(cubeSize, 128) == 0));

assert(all(mod(bbox(:, 1), 128) == 1));

X = [bbox(1, 1):cubeSize(1):bbox(1, 2), bbox(1, 2) + 1];
Y = [bbox(2, 1):cubeSize(2):bbox(2, 2), bbox(2, 2) + 1];
Z = [bbox(3, 1):cubeSize(3):bbox(3, 2), bbox(3, 2) + 1];

dimCount = cellfun(@numel, {X, Y, Z}) - 1;
inputCell = cell(prod(dimCount), 1);

for i=1:dimCount(1)
    for j=1:dimCount(2)
        for k=1:dimCount(3)
            idx = sub2ind( ...
                [length(X)-1 length(Y)-1 length(Z)-1], i, j, k);
            inputCell{idx} = {[X(i) X(i+1)-1; Y(j) Y(j+1)-1; Z(k) Z(k+1)-1]};
        end
    end
end
    
    % init wkw dataset, if needed
if isfield(p.class, 'backend') && strcmp(p.class.backend, 'wkwrap')
	wkwInit('new', p.classMyelin.root, 32, 32, 'single', 1);
end

if p.cnn.GPU
    numGpus = 1;
else
    numGpus = 0;
end

job = Cluster.startJob( ...
    @jobWrapper, inputCell, ...
    'sharedInputs', {p}, ...
    'name', 'myelinDetection', ...
    'cluster', { ...
        'memory', 12, ...
        'time', '24:00:00', ...
        'gpus', numGpus});
end

function jobWrapper(p, bbox)
    onlyFwdPass3DonKnossosFolder( ...
        p.cnn.myelin, p.cnn.GPU, p.raw, ...
        p.classMyelin, bbox, p.norm.func);
end

