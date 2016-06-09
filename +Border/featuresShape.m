function [weights, weightNames] = featuresShape(params, voxelIds)
    % Written by
    %   Benjamin Ehret <no-email-known>
    %   Manuel Berning <manuel.berning@brain.mpg.de>
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    
    voxelIdsSize = size(voxelIds);
    voxelIds = voxelIds(:);
    
    % prepare 
    cubeSize = params.tileSize;
    doFunc = @(ids) main(cubeSize, ids);
    
    % do it!
    weights = cellfun( ...
        doFunc, voxelIds, ...
        'UniformOutput', false);
    
    % reshape output
    weights = reshape( ...
        weights, voxelIdsSize);
    
    weightNames = {
        'shape: logObjSize', ...
        'shape: pcFracs(1)', ...
        'shape: pcFracs(2)', ...
        'shape: pcFracs(3)', ...
        'shape: pcFracs(3)/pcFracs(2)', ...
        'shape: pcFracs(3)/pcFracs(1)', ...
        'shape: pcFracs(2)/pcFracs(1)', ...
        'shape: entropy' };
end

function weights = main(cubeSize, voxelIds)
    voxelCount = numel(voxelIds);
    voxelList = nan(voxelCount, 3);

    [voxelList(:, 1), voxelList(:, 2), voxelList(:, 3)] = ...
        ind2sub(cubeSize, voxelIds);
    
    % Hmm not sure why Benjamin chose PC (and only latent output), extend?
    [~, ~, latent] = pca( ...
        voxelList, 'Economy', true);
    pcFracs = latent ./ sum(latent);
    
    % prepare output
    weights = zeros(1, 8);
    
    % skip small objects
    if objSize <= 3; return; end;
    
    % calculate features
    weights(1) = log(objSize);
    weights(2) = pcFracs(1);
    weights(3) = pcFracs(2);
    weights(4) = pcFracs(3);
    
    weights(5) = pcFracs(3) / pcFracs(2);
    weights(6) = pcFracs(3) / pcFracs(1);
    weights(7) = pcFracs(2) / pcFracs(1);
    
    % entropy
    weights(8) = ...
        -pcFracs(1) * log(pcFracs(1)) ...
        -pcFracs(2) * log(pcFracs(2)) ...
        -pcFracs(3) * log(pcFracs(3));
    
    % replate NaN values by zero
    weights(isnan(weights)) = 0;
end