function [featVals, featNames] = ...
        featuresTexture(param, box, voxelIds)
    % featuresTexture(param, cubeIdx, voxelIds)
    %   Calculates the texture features for all the regions
    %   specified in 'voxelIds'.
    %
    % param
    %   Parameters produced by 'setParameterSettings.m'
    %
    % box
    %   Bounding box (in global coordinates) around the region
    %   for which the filter responses should be calculated.
    %
    % voxelIds
    %   Cell array. Each entry contains the list of linear
    %   indices of all voxels making up the region.
    %
    % Written by
    %   Manuel Berning <manuel.berning@brain.mpg.de>
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    % config
    padSize = 10;

    % add padding to avoid border effects
    box(:, 1) = box(:, 1) - padSize;
    box(:, 2) = box(:, 2) + padSize;
    
    % reshape borders
    voxelIdsSize = size(voxelIds);
    voxelIds = voxelIds(:);

    % prepare output
    featVals = nan(numel(voxelIds), 0);
    featNames = cell(0);

    inputs = param.feature.input;
    filters = param.filter;

    for curInputIdx = 1:numel(inputs)
        % load input data
        curInput = param.feature.input{curInputIdx};
        data = loadInputData(param, curInput, box);

        for curFilterIdx = 1:numel(filters)
            curFilter = filters{curFilterIdx};

            % filter
            curName = curFilter{1};
            curSizes = curFilter{2};
            curMore = curFilter{3};

            for curSizeIdx = 1:numel(curSizes)
                curSize = curSizes(curSizeIdx);

                % apply filter
                filteredData = filter3d( ...
                    curName, data, curSize, curMore);

                % make cell array
                if ~iscell(filteredData)
                    filteredData = {filteredData};
                end

                for p = 1:numel(filteredData)
                    % remove padding
                    filteredData{p} = filteredData{p}( ...
                        (padSize + 1):(end - padSize), ...
                        (padSize + 1):(end - padSize), ...
                        (padSize + 1):(end - padSize));

                    % sanity check
                    assert(isreal(filteredData{p}));

                    % compute features
                    [newFeatVals, newFeatNames] = ...
                        featureDesign(filteredData{p}, voxelIds);

                    % build feature names
                    newFeatNamePrefix = sprintf( ...
                        '%s %s filtersize=%d feature=%d ', ...
                        curInput, curName, curSize, curMore);
                    newFeatNames = cellfun(@(str) ...
                        {[newFeatNamePrefix, str]}, newFeatNames);

                    % add features to output
                    featVals = [featVals, newFeatVals];
                    featNames = [featNames(:)', newFeatNames(:)'];
                end
            end
        end
    end
    
    % restore shape of borders
    featVals = num2cell(featVals, 2);
    featVals = reshape(featVals, voxelIdsSize);
end

function data = loadInputData(param, inputName, box)
    switch inputName
        case 'raw'
            % raw data
            data = loadRawData( ...
                param.raw.root, param.raw.prefix, box);
            data = param.norm.func(single(data));
            
        case 'aff'
            % affinity / membrane
            data = loadClassData( ...
                param.class.root, param.class.prefix, box);
            
        case 'mito'
            % mitochondria
            
            % TODO
            %   factor out into parameter structure
            mitoDataRoot = ...
                '/gaba/u/mberning/results/pipeline/20141007T094904/mvm/';
            mitoDataPrefix = ...
                '2012-09-28_ex145_07x2_corrected_mag1';
            data = loadMitoData( ...
                mitoDataRoot, mitoDataPrefix, box);
            
        otherwise
            error('Invalid feature input');
            
    end
end