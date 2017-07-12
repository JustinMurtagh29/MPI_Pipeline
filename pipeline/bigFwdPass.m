function job = bigFwdPass( p, bbox )
    % job = bigFwdPass( p, bbox )
    %   Applies the CNN for membrane detection to the raw data
    %   within the specified bounding box.
    
    % sanity check
    assert(all(mod(bbox(:, 1), 128) == 1));
    
    % This is not of any importance due to CNN translation invariance,
    % can be choosen for computational efficency, currenlty optimized for
    % running on GPU with 12GB, should be multiples of 128, this is same
    % as tileSize right now, no reason it has to be.
    cubeSize = [512 512 256];
    assert(all(mod(cubeSize, 128) == 0));
    
    X = [bbox(1, 1):cubeSize(1):bbox(1, 2), p.bbox(1, 2) + 1];
    Y = [bbox(2, 1):cubeSize(2):bbox(2, 2), p.bbox(2, 2) + 1];
    Z = [bbox(3, 1):cubeSize(3):bbox(3, 2), p.bbox(3, 2) + 1];
    
    dimCount = cellfun(@numel, {X, Y, Z}) - 1;
    inputCell = cell(prod(dimCount), 1);
    
    for i=1:dimCount(1)
        for j=1:dimCount(2)
            for k=1:dimCount(3)
                idx = sub2ind( ...
                    [length(X)-1 length(Y)-1 length(Z)-1], i, j, k);
                inputCell{idx} = { ...
                    p.cnn.first, p.cnn.GPU, p.raw, p.class, ...
                    [X(i) X(i+1)-1; Y(j) Y(j+1)-1; Z(k) Z(k+1)-1], p.norm.func};
            end
        end
    end
    
    % init wkw dataset, if needed
    if isfield(p.class, 'backend') && strcmp(p.class.backend,'wkwrap')
        wkwInit('new', p.class.root, 32, 32, 'single', 1);
    end

    % Same function for all input arguments
    functionH = @onlyFwdPass3DonKnossosFolder;

    if p.cnn.GPU
        job = startGPU(functionH, inputCell, 'classification');
    else
        job = startCPU(functionH, inputCell, 'classification');
    end
end

