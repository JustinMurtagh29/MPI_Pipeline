function [weights, weightNames] = featureDesign(input, voxelIds)
    % [weights, weightNames] = featureDesign(input, voxelIds)
    %   Computes weights / features / statistics for each
    %   region defined by 'voxelIds'.
    %
    % input
    %   NxMxL matrix of real numbers. Filter response.
    %
    % voxelIds
    %   Cell array. Each entry contains the linear indices
    %   of all voxels making up a region.
    %
    % Written by
    %   Manuel Berning <manuel.berning@brain.mpg.de>
    
    % sanity check
    assert(isreal(input));

    % Finding value for each voxel
    voxelVals = cellfun( ...
        @(ids) input(ids(:))', voxelIds(:), ...
        'UniformOutput', false);

    % Computing 7 features
    % - 5 Quantiles: 0, 0.25, 0.5, 0.75, 1
    weights = cell2mat(cellfun( ...
        @(t) quantile(t, [0, 0.25, 0.5, 0.75, 1]), ...
        voxelVals, 'UniformOutput', false));

    % - Mean
    weights = [weights, cellfun(@mean, voxelVals)];

    % - Standard deviation
    weights = [weights, cellfun(@std, voxelVals)];

    weightNames = { ...
        'quantile0', ...
        'quantile0.25', ...
        'quantile0.5', ...
        'quantile0.75', ...
        'quantile1', ...
        'mean', ...
        'std'};
end

