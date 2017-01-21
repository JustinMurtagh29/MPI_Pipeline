function job = globalizeCorrespondences(p)
    % Globalize correspondences
    files = dir([p.correspondence.saveFolder, '*.mat']);

    % Bad but internal: Sort out global correspondences in case pipeline has already run
    % amotta: why the == 16? what does it mean?
    idx = cellfun(@length, {files(:).name}) == 16;
    files = files(idx);

    % collect parameters
    fileCount = numel(files);
    inputCell = cell(fileCount, 1);

    for i = 1:fileCount
        inputCell{i} = {p, files(i).name};
    end

    functionH = @globalCorrSeg;
    job = startCPU(functionH, inputCell, 'globalCorrespondences', 12, 100);
end
