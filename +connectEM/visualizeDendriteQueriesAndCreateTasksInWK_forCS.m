% Load p structure for bounding box of the dataset and transfer to wK coordinates
    load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat');
    
    % Where to find query information (position and euler angles)
    batchFolder = '/gaba/u/mberning/results/pipeline/20170217_ROI/aggloState/dendriteQueries_3/';

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
    
    % Set bounding box area around query position
    shift10um = round([10*1000 10*1000 10*1000]./[11.24 11.24 28]);

    for i=1:100
        batch = load(fullfile(batchFolder, ['batch' num2str(i, '%.4i') '.mat']));
        for j=1:numel(batch.q.pos)
            for k=1:numel(batch.q.pos{j})
                % Calculate bbox
                bbox(1:3,1) = batch.q.pos{j}{k}(:) - shift10um';
                bbox(4:6,1) = batch.q.pos{j}{k}(:) + 2*shift10um';
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

