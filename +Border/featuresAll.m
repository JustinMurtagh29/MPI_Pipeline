function featuresAll(param, cubeIdx)
    % featuresAll(param, cubeIdx)
    %   Calculates the features (texture and shape thus far)
    %   for all extended borders in a specified segmentation
    %   cube.
    %
    % param
    %   Parameters returned by 'setParameterSettings'
    %
    % cubeIdx
    %   Linear index of current cube
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    % get important parameters
    cubeParams = param.local(cubeIdx);
    cubeDir = cubeParams.saveFolder;

    % load pixel list for borders and segments
    load([cubeDir, 'bordersExt.mat'], 'borders');
    
    disp('Shape features...');
    [shapeFeatVals, shapeFeatNames] = ...
        Border.featuresShape(param.tileSize, borders);
    
    disp('Texture features...');
    [texFeatVals, texFeatNames] = ...
        Border.featuresTexture(param, cubeIdx, borders);
    
    % sanity check
    assert(all(size(texFeatVals) == size(shapeFeatVals)));
    
    % build output
    featVals = cellfun( ...
        @(tVals, sVals) [tVals(:)', sVals(:)'], ...
        texFeatVals, shapeFeatVals, 'UniformOutput', false);
    featNames = [texFeatNames(:)', shapeFeatNames(:)'];
    
    % prepare output
    outStruct = struct;
    outStruct.feats = featVals;
    outStruct.featNames = featNames;
    
    outFile = [cubeDir, 'bordersExtFeats.mat'];
    save(outFile, '-struct', 'outStruct');
end