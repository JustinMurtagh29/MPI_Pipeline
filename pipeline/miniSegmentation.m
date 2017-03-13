function job = miniSegmentation(p)

% allocate task inputs
taskCount = numel(p.local);
inputCell = cell(taskCount, 1);

for idx = 1:taskCount
    % create cube folder, if needed
    if ~exist(p.local(idx).saveFolder, 'dir')
        mkdir(p.local(idx).saveFolder);
    end
    
    if isfield(p.local(idx), 'class')
        % This is for pT train, probably needs update at some point?
        inputCell{idx} = { ...
            p.local(idx).class.root, p.local(idx).class.prefix, ...
            p.local(idx).bboxBig, p.seg.func, p.local(idx).segFile};
    else
        inputCell{idx} = { ...
            p.class.root p.class.prefix, ....
            p.local(idx).bboxBig, p.seg.func, p.local(idx).tempSegFile};
    end
end

functionH = @segmentForPipeline;
job = startCPU(functionH, inputCell, 'segmentation', 36);

end

