function job = miniSegmentation(p)
inputCellBuilder = @(idx) {{ ...
    p.class, p.local(idx).bboxBig, ...
    p.seg.func, p.local(idx).tempSegFile}};

inputCell = arrayfun(inputCellBuilder, 1:numel(p.local));
job = startCPU(@segmentForPipeline, inputCell, 'segmentation', 36);
end
