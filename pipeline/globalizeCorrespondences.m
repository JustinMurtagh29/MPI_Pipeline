function job = globalizeCorrespondences(p)
% Globalize correspondences
    
    files = dir([p.correspondence.saveFolder, '*.mat']);
    % Bad but internal: Sort out global correspondences in case pipeline has already run
    idx = cellfun(@length, {files(:).name}) == 16;
    files = files(idx);
    for i=1:length(files)
        inputCell{i} = {p, files(i).name};
    end
    functionH = @globalCorrSeg;
    job = startCPU(functionH, inputCell, 'globalCorrespondences', 12, 100);



end
