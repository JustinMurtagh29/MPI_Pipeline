function job = miniSegmentation(p)

% allocate task inputs
taskCount = numel(p.local);
inputCell = cell(taskCount, 1);

for idx = 1:taskCount
    % create cube folder, if needed
    if ~exist(p.local(idx).saveFolder, 'dir')
        mkdir(p.local(idx).saveFolder);
    end
    
    inputCell{idx} = { ...
        p.class, p.local(idx).bboxBig, ...
        p.seg.func, p.local(idx).tempSegFile};
end

functionH = @segmentForPipeline;
job = startCPU(functionH, inputCell, 'segmentation', 36);

end

