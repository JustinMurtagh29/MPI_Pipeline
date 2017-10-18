% Load p structure for bounding box of the dataset and transfer to wK coordinates
    load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat');
    % Switch to wk-coordinates (-1) and extend border by 128 so that CNN border region
    % is also shown (will lead to slightly larger bounding box
    bbox = bsxfun(@minus, p.bbox , [1; 1; 1]) + repmat([-128 128],3,1);
    % Switch to width height depth format
    bbox(:,2) = bbox(:,2) - bbox(:,1) + 1;
    % Linearize
    bbox = reshape(bbox, 6, 1);

    % Where to find query information (position and euler angles)
    batchFolder = '/gaba/u/mberning/results/pipeline/20170217_ROI/aggloState/dendriteQueries_2/';

    % Define constant task parameters
    dataSet = '2012-09-28_ex145_07x2_20171017_DendQueries';
    taskTypeId = '5936fb14760000e30a906635';
    experienceDomain = 'focus_flight';
    minExperience = 1;
    % Positions and rotations loaded in script below
    instances = 1;
    team = 'Tracing crew';
    project = 'CS_MB_dendriteEndingQueries_20171017';
    scriptId = '591c8d496e000012577043a3';

    % Where to write results
    outputFile = fullfile(batchFolder, 'queries.txt');
    fileID = fopen(outputFile,'w');

    for i=1:100
        batch = load(fullfile(batchFolder, ['batch' num2str(i, '%.4i') '.mat']));
        for j=1:numel(batch.q.pos)
            for k=1:numel(batch.q.pos{j})
                taskString = sprintf('%s %s %s %u %u %u %u %.2f %.2f %.2f %u %s %u %u %u %u %u %u %s %s', ...
                    dataSet, ...
                    taskTypeId, ...
                    experienceDomain, ...
                    minExperience, ...
                    batch.q.pos{j}{k}(:), ...
                    batch.q.angles{j}{k}(:), ...
                    instances, ...
                    team, ...
                    bbox(:), ...
                    project, ...
                    scriptId);
                fprintf(fileID, '%s\n', taskString);
            end
        end
    end

    fclose(fileID);

