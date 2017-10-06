function output = generateQueriesFromChiasmata( ...
        aggloFile, chiasmataId, outputFolder)
    % load axon agglomerates
    agglos = load(aggloFile, 'axons', 'indBigAxons');
    
    curDateStr = datestr(clock, 30);
    fid = sprintf('%s_flightTasks.txt', curDateStr);
    fid = fopen(fullfile(outputFolder, fid), 'w');
    
    output = [];
    for agglo_idx = find(agglos.indBigAxons)'
        temp = load(fullfile( ...
            '/tmpscratch/kboerg/chiasmata', ...
            sprintf('chiasmataX%s_%d', chiasmataId, floor(agglo_idx / 100)), ...
            sprintf('visX%s_%d', chiasmataId, agglo_idx), ...
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
                extend = round(4000 ./[11.24,11.24,28]);
                dd = sqrt(sum((temp.output.direction{i}(j,:) .* [11.24,11.24,28]).^2));
                if dd>2000
                    extend = round(2*dd ./[11.24,11.24,28]);
                end
                [phi, theta, psi] = calculateEulerAngles(-temp.output.direction{i}(j,:), [11.24,11.24,28]);
                q.angles{i}(j,:) = [phi theta psi];
                minPos = temp.output.position{i}(j,:) - extend;
                sizeBbox = 2*extend;
                taskString = ['2012-09-28_ex145_07x2_ROI2017,56d6a7c6140000d81030701e,queriesMHchiasma,1,' ...
                    num2str(temp.output.position{i}(j,1)-1) ',' num2str(temp.output.position{i}(j,2)-1) ',' num2str(temp.output.position{i}(j,3)-1) ',' ...
                    num2str(phi) ',' num2str(theta) ',' num2str(psi) ',1,Connectomics department,' ...
                    num2str(minPos(1)) ',' num2str(minPos(2)) ',' num2str(minPos(3)) ',' ...
                    num2str(sizeBbox(1)) ',' num2str(sizeBbox(2)) ',' num2str(sizeBbox(3)) ',' 'queriesMHchiasma'];
                fprintf(fid, '%s\n', taskString);
                output(end+1,:) = [agglo_idx, i, j, temp.output.position{i}(j,:), temp.output.ccCenterIdx(i)];
            end
        end
    end
    
    fclose(fid);
    
    % NOTE(amotta): Store `output` on disk so that we can avoid running all
    % of the above stuff in the future (e.g., when splitting chiasmata).
    out = struct;
    out.aggloFile = aggloFile;
    out.chiasmataId = chiasmataId;
    
    out.queries = output;
    out.gitInfo = Util.gitInfo();
    
    saveFile = sprintf('%s_data.mat', curDateStr);
    saveFile = fullfile(outputFolder, saveFile);
    Util.saveStruct(saveFile, out);
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
