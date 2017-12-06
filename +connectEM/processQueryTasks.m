function processQueryTasks(i, batchBoundaries, axons, borderPositions, ...
    directions, p, outputFolder)

    % Written by
    %   Manuel Berning <manuel.berning@brain.mpg.de>
    %   Christian Schramm <christian.schramm@brain.mpg.de>
    % Collects the information for querying agglo endings. The number of 
    % agglos is presumly (connectEM.generateAxonQueries) divided into 100
    % batches for parallelisation.
    % Input parameters:
    %   i: int
    %     index of current batch
    %   batchBoundaries
    %     (1xN) array containing the start indices for each batch. The last
    %       index for a batch is given by batchBoundaries(i+1)-1, i.e.
    %       batchBoundaries(end+1) should be length(axons) + 1
    %   axons
    %     axon agglomerates in the old format (i.e. list of segment ids)
    %   borderPositions
    %     border positions of chosen endings
    %   directions
    %     pca of border ID
    %   p: struct
    %     segmentation parameter file
    %   outputFolder: string
    %     output directory. output files are of the format 'batch_\d.mat'

    theseIdx = batchBoundaries(i):batchBoundaries(i+1)-1;
    theseAxons = axons(theseIdx);
    theseBorderPositions = borderPositions(theseIdx);
    theseDirections = directions(theseIdx);

    q.pos = {};
    q.posBorder = {};
    q.dir = {};
    q.angles = {};
    for j=1:length(theseAxons)
        tic;
        q.posBorder{end+1} = theseBorderPositions{j};
        thesePositions = cellfun(@(x,y,z) ...
            connectEM.correctQueryLocationToEndOfSegment(p, x, y, z, 200), ...
            repmat(theseAxons(j),size(theseBorderPositions{j},1),1), ...
            theseBorderPositions{j}, ...
            theseDirections{j}, 'uni', 0);
        q.pos{end+1} = thesePositions;
        q.dir{end+1} = theseDirections{j};
        [phi, theta, psi] = cellfun(@(x) ...
            connectEM.calculateEulerAngles(x, p.raw.voxelSize), ...
            q.dir{end});
        q.angles{end+1} = mat2cell(cat(2, phi, theta, psi), ...
            ones(numel(phi),1), 3);
        toc;
    end

    % Calculate euler angles & put into old format
    save(fullfile(outputFolder, ['batch' num2str(i, '%.4i') '.mat']), ...
        'q', 'theseAxons');
    display(['Batch ' num2str(i, '%.4i') ' done']);
    clear these* q phi theta psi;
end

