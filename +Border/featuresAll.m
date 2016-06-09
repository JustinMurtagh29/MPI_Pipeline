function featuresAll(param, cubeIdx)
    % get important parameters
    cubeParams = param.local(cubeIdx);
    cubeDir = cubeParams.saveFolder;

    % load pixel list for borders and segments
    load([cubeDir, 'bordersExt.mat'], 'borders');
    
    disp('Shape features...');
    [shapeFeatVals, shapeFeatNames] = ...
        Border.featuresShape(param, borders);
    
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