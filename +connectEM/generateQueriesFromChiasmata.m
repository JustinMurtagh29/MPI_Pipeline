function output = generateQueriesFromChiasmata( ...
        aggloFile, chiasmataId, outputFolder)
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
            '/tmpscratch/kboerg/chiasmata', ...
            sprintf('chiasmataX%s_%d', chiasmataId, floor(aggloIdx / 100)), ...
            sprintf('visX%s_%d', chiasmataId, aggloIdx), ...
            'result.mat'));
        
        for i = 1:numel(temp.output.position)
            % sanity check
            assert(size(temp.output.direction{i}, 1) ...
                == temp.output.nrExits(temp.output.ccCenterIdx(i)));
            curNrExits = size(temp.output.direction{i}, 1);
            
            if curNrExits < 4
                curNrQueries = 0;
            elseif curNrExits == 4
                curNrQueries = 1;
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
                    temp.output.ccCenterIdx(i)]; %#ok
                taskDef(end + 1, {'position', 'rotation', 'bbox'}) = { ...
                    temp.output.position{i}(j, :) - 1, ...
                   [phi, theta, psi], [minPos, sizeBbox]}; %#ok
            end
        end
    end
    
    % NOTE(amotta): Store `output` on disk so that we can avoid running all
    % of the above stuff in the future (e.g., when splitting chiasmata).
    out = struct;
    out.aggloFile = aggloFile;
    out.chiasmataId = chiasmataId;
    
    out.queries = output;
    out.gitInfo = Util.gitInfo();
    
    saveFile = sprintf('%s_data.mat', curDateStr);
    Util.saveStruct(fullfile(outputFolder, saveFile), out);
    
    taskDefFile = sprintf('%s_flightTasks.txt', curDateStr);
    writeTaskDefinition(taskDef, fullfile(outputFolder, taskDefFile));
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

function writeTaskDefinition(taskDef, outFile)
    fullTaskDef = taskDef;
    fullTaskDef.dataSet(:) = {'2012-09-28_ex145_07x2_ROI2017'};
    fullTaskDef.taskTypeId(:) = {'56d6a7c6140000d81030701e'};
    fullTaskDef.expDomain(:) = {'queriesMHchiasma'};
    fullTaskDef.expMinVal(:) = 1;
    fullTaskDef.instances(:) = 1;
    fullTaskDef.team(:) = {'Connectomics department'};
    fullTaskDef.project(:) = {'queriesMHchiasma'};
    
    % reorder-columns
    columnOrder = { ...
        'dataSet', 'taskTypeId', 'expDomain', 'expMinVal', ...
        'position', 'rotation', 'instances', 'team', 'bbox', 'project'};
    fullTaskDef = fullTaskDef(:, columnOrder);
    
    % write CSV file
    writetable(fullTaskDef, outFile, 'WriteVariableNames', false);
end
