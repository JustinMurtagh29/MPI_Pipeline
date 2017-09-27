function output = generateQueriesFromChiasmata(outputFolder, temp2)
    % TODO change hardcoded paths and ==4 vs >4 distinction
    fid = fopen([outputFolder datestr(clock,30) '_flightTasks.txt'], 'w');
    output = [];
    for idx_agglo = find(temp2.indBigAxons)'
        idx_agglo
        temp = load(['/tmpscratch/kboerg/chiasmata/chiasmataX33_' num2str(floor(idx_agglo/100)) '/visX33_' num2str(idx_agglo) '/result.mat']);
        for i=1:length(temp.output.position)
            if temp.output.nrExits(temp.output.ccCenterIdx(i))<4 
                continue
            end
            for j=1 % 2:size(temp.output.direction{i},1) %only one task per chiasma CHANGE FOR >4 CHIASMATA
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
                output(end+1,:) = [idx_agglo, i, j, temp.output.position{i}(j,:), temp.output.ccCenterIdx(i)];
            end
        end
    end
    fclose(fid);
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
