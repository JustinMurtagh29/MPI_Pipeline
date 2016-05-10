function collectGlobalCoMs(p, maxID)
    % Calculate & collect all CoMs of segments in supervoxel graph

    % Submit job to jobmanager for calculation
    idx = 1;
    for i=1:size(p.local,1)
        for j=1:size(p.local,2)
            for k=1:size(p.local,3)
                inputCell{idx} = {p, i, j, k}; 
                idx = idx + 1;
            end
        end
    end
    functionH = @calculateLocalCoMs;
    job = startCPU(functionH, inputCell, 'CoM calculation in small bounding box');

    wait(job, 'finished');

    % Collect
    cubesProcessed = 0;
    com = zeros(maxID, 3);
    tic;
    for i=1:size(p.local,1)
        for j=1:size(p.local,2)
            for k=1:size(p.local,3)
                load([p.local(i,j,k).saveFolder 'comSmall.mat']);
                com(minObjId:(minObjId+lengthRP-1),:) = globalPos;
                % Display progess
                cubesProcessed = cubesProcessed + 1;
                Util.progressBar(cubesProcessed, prod(size(p.local)));
            end
        end
    end
    save([p.saveFolder 'CoM.mat'], 'com');

end

