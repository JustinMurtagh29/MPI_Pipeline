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

functionH = @miniSegmentationJobWrapper;
job = startCPU(functionH, inputCell, 'segmentation', 36);

end

function miniSegmentationJobWrapper(saveFolder, idx)
% INPUT saveFolder: string
%           Path to segmentation main folder containing the
%           'allParameter.mat' file.
%       idx: int
%           Linear index of the current local cube.

m = load(fullfile(saveFolder, 'allParameter.mat'));
p = m.p;
segmentForPipeline(p.class.root, p.class.prefix, p.local(idx).bboxBig, ...
    p.seg.func, p.local(idx).tempSegFile);
end

