function output = generateQueriesFromChiasmata( ...
        aggloFile, chiasmaDir, chiasmaId, outputFolder)
    % Written by
    %   Kevin Boergens <kevin.boergens@brain.mpg.de>
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    % TODO(amotta): Use parameter structure
    voxelSize = [11.24, 11.24, 28];
    curDateStr = datestr(clock, 30);
    
    % load axon agglomerates
    agglos = load(aggloFile, 'axons', 'indBigAxons');
    
    output = [];
    taskDef = table;
    taskDef.position = zeros(0, 3);
    taskDef.rotation = zeros(0, 3);
    taskDef.bbox = zeros(0, 6);
    
    for aggloIdx = find(agglos.indBigAxons)'
        temp = load(fullfile( ...
            chiasmaDir, ...
            sprintf('chiasmataX%s_%d', chiasmaId, floor(aggloIdx / 100)), ...
            sprintf('visX%s_%d', chiasmaId, aggloIdx), ...
            'result.mat'));
        
        for i = 1:numel(temp.output.position)
            % sanity check
            assert(size(temp.output.direction{i}, 1) ...
                == temp.output.nrExits(temp.output.ccCenterIdx(i)));
            curNrExits = size(temp.output.direction{i}, 1);
            
            if curNrExits < 4
                curNrQueries = 0;
            else
                curNrQueries = curNrExits;
            end
            
            for j = 1:curNrQueries
                extend = round(4000 ./ voxelSize);
                dd = sqrt(sum((voxelSize ...
                    .* temp.output.direction{i}(j,:)) .^ 2));
                if dd > 2000
                    extend = round(2 * dd ./ voxelSize);
                end
                
                % determine size of bounding box
                minPos = temp.output.position{i}(j, :) - extend;
                sizeBbox = 2*extend;
                
                % determine Euler angles
                [phi, theta, psi] = calculateEulerAngles( ...
                    -temp.output.direction{i}(j, :), voxelSize);
                output(end+1, :) = [ ...
                    aggloIdx, i, j, ...
                    temp.output.position{i}(j, :), ...
                    temp.output.ccCenterIdx(i), ...
                    temp.output.queryIdx{i}(j)]; %#ok
                taskDef(end + 1, {'position', 'rotation', 'bbox'}) = { ...
                    temp.output.position{i}(j, :) - 1, ...
                   [phi, theta, psi], [minPos, sizeBbox]}; %#ok
            end
        end
    end

    rng(0);
    randIds = randperm(size(output, 1));
    
    % shuffle chiasmata
    output = output(randIds, :);
    taskDef = taskDef(randIds, :);
    
    % make batches of 500 chiasmata
   [~, ~, uniChiasma] = unique(output(:, 1:2), 'rows', 'stable');
    
    for curLimit = 1:500:max(uniChiasma)
        curMask = uniChiasma < curLimit;
        
        % bubble to the top
        output = cat(1, output(curMask, :), output(~curMask, :));
        taskDef = cat(1, taskDef(curMask, :), taskDef(~curMask, :));
        uniChiasma = cat(1, uniChiasma(curMask), uniChiasma(~curMask));
    end
    
    % NOTE(amotta): Store `output` on disk so that we can avoid running all
    % of the above stuff in the future (e.g., when splitting chiasmata).
    out = struct;
    out.aggloFile = aggloFile;
    out.chiasmataId = chiasmaId;
    
    out.queries = output;
    out.taskDef = taskDef;
    out.gitInfo = Util.gitInfo();
    
    saveFile = sprintf('%s_data.mat', curDateStr);
    Util.saveStruct(fullfile(outputFolder, saveFile), out);
    
    taskParam = struct;
    taskParam.dataSet = '2012-09-28_ex145_07x2_ROI2017';
    taskParam.taskTypeId = '56d6a7c6140000d81030701e';
    taskParam.expDomain = 'queriesMHchiasma';
    taskParam.expMinVal = 1;
    taskParam.instances = 1;
    taskParam.team = 'Connectomics department';
    taskParam.project = 'queriesMHchiasma';
    
    taskDefFile = sprintf('%s_flightTasks.txt', curDateStr);
    connectEM.exportTaskDefinitions( ...
        taskParam, taskDef, fullfile(outputFolder, taskDefFile));
end

function [phi, thetha, psi] = calculateEulerAngles(di, voxelSize)
    % Calculate angles (as deinfed in wK) in degrees from direction vector

    % Make sure direction is normalized
    di = di .* voxelSize;
    di = di ./ norm(di);
    % Calculate euler angles, from wk-paper -> RESCOPaaS
    eulerAngles = connectEM.diffToEulerAngle(di);
    phi = round(eulerAngles(1));
    thetha = round(eulerAngles(2));
    psi = round(eulerAngles(3));
end
