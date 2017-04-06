function job = globalizeCorrespondences(p)
    % Globalize correspondences
    files = dir([p.correspondence.saveFolder, '*.mat']);

    % Bad but internal: Sort out global correspondences in case pipeline has already run
    % amotta: why the == 16? what does it mean?
    % Length of the file name, here we only want the mat files without global suffix in case this function ran before
    % Agreed that this is suboptimal
    idx = cellfun(@length, {files(:).name}) == 16;
    files = files(idx);

    % collect parameters
    fileCount = numel(files);
    sharedInputs = {p};
    inputCell = cell(fileCount, 1);

    for i = 1:fileCount
        inputCell{i} = {files(i).name};
    end

    functionH = @globalCorrSeg;
    job = Cluster.startJob( ...
        functionH, inputCell, ...
        'name', 'globalCorrespondences', ...
        'sharedInputs', sharedInputs, ...
        'cluster', '-l h_vmem=12G', ...
        'taskGroupSize', 100);
end
