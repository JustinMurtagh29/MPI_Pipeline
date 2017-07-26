function raw = copyDebugRoi()
    % Copy part of 2012_09_28_ex145_07x2_new2 dataset availible on webKnossos
    % to a new dataset to have the ROI for segmentation aligned to KNOSSOS
    % cubes. Also, the copied ROI is such that we generate partial cubes at
    % the end. This should allow us to debug the pipeline.
    
    % config
    center = 1 + [3811, 7860, 1866];
    roiSize = [512, 512, 256] + [97, 97, 97];
    cnetBorder = [25, 25, 10];
    
    % Define source of copy process
    src = struct;
    src.root = '/gaba/u/kboerg/st07x2_new2/color/1/';
    src.prefix = '2012-09-28_ex145_07x2_new2_mag1';
    
    srcBox = nan(3, 2);
    srcBox(:, 1) = round(center - roiSize / 2);
    srcBox(:, 2) = srcBox(:, 1) + roiSize(:) - 1;
    
    srcBox(:, 1) = srcBox(:, 1) - cnetBorder(:);
    srcBox(:, 2) = srcBox(:, 2) + cnetBorder(:);
    
    % Define destination of copy process
    dest = struct;
    dest.root = '/u/amotta/pipeline/2012-09-28_ex145_07x2_ROI2016/data/knossos/raw/';
    dest.prefix = '2012-09-28_ex145_07x2_ROI2016_mag1';
    destOffset = 129 - cnetBorder;

    % Read from one, write to wKcubes
    raw = loadRawData(src, srcBox);
    saveRawData(dest, destOffset, raw);
end
