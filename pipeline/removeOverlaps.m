function removeOverlaps(p)

    % Transfer local cubes to segFile
    for i=1:numel(p.local)
        inputCell{i} = {p.local(i).tempSegFile, p.local(i).segFile, p.local(i).bboxSmall, p.local(i).bboxBig};
    end
    functionH = @rewriteWithoutOverlaps;
    job = startCPU(functionH, inputCell, 'overlapRemoval', 12, 10);
    Cluster.waitForJob(job);

    % Modify all correspondences accordingly
    display('Modifying correspondences');
    saveFolder = p.correspondence.saveFolder;
    files = dir([saveFolder '*.mat']);
    tic;
    for i=1:length(files)
        % Exclude *global.mat written in a later step if already present (reruns)
        if isempty(strfind(files(i).name, 'global.mat'))
            % Load correspondences
            corr = load([saveFolder files(i).name]);
            result = renumberCorrespondences(p, corr);
            Util.saveStruct([saveFolder files(i).name], result); 
        end
        Util.progressBar(i, length(files));
    end
end

function new = renumberCorrespondences(p, old)
    % Load segment renumbering
    idx = old.cubeCoords1;
    cube1 = load(p.local(idx(1),idx(2),idx(3)).segFile, 'oldSegments', 'newSegments');
    idx = old.cubeCoords2;
    cube2 = load(p.local(idx(1),idx(2),idx(3)).segFile, 'oldSegments', 'newSegments');
    new = old;
    new.uniqueCorrespondences = modifyCorrespondences(new.uniqueCorrespondences, cube1, cube2);
    new.doubleCorrespondences = modifyCorrespondences(new.doubleCorrespondences, cube1, cube2);
end

function newCorr = modifyCorrespondences(corr, cube1, cube2)
    newCorr(:,1) = changem(corr(:,1), cube1.newSegments, cube1.oldSegments);
    newCorr(:,2) = changem(corr(:,2), cube2.newSegments, cube2.oldSegments);
end

