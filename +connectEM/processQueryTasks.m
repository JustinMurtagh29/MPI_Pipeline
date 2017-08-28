function processQueryTasks(i,batchBoundaries,axons,borderPositions,directions,p,outputFolder)

    theseIdx = batchBoundaries(i):batchBoundaries(i+1)-1;
    theseAxons = axons(theseIdx);
    theseBorderPositions = borderPositions(theseIdx);
    theseDirections = directions(theseIdx);

    q.pos = {};
    q.dir = {};
    q.angles = {};
    for j=1:length(theseAxons)
        tic;
        thesePositions = cellfun(@(x,y,z)connectEM.correctQueryLocationToEndOfSegment(p, x, y, z, 200), ...
            repmat(theseAxons(j),size(theseBorderPositions{j},1),1), theseBorderPositions{j}, theseDirections{j}, 'uni', 0);
        q.pos{end+1} = thesePositions;
        q.dir{end+1} = theseDirections{j}(~cellfun('isempty',thesePositions));
        [phi, theta, psi] = cellfun(@(x)connectEM.calculateEulerAngles(x, p.raw.voxelSize), q.dir{end});
        q.angles{end+1} = mat2cell(cat(2, phi, theta, psi), ones(numel(phi),1), 3);
        toc;
    end

    % Calculate euler angles & put into old format
    save([outputFolder 'batch' num2str(i, '%.4i') '.mat'], 'q', 'theseAxons');
    display(['Batch ' num2str(i, '%.4i') ' done']);
    clear these* q phi theta psi;
        