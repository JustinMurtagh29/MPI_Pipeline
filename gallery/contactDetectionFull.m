function contactDetectionFull()

outputDir = '/zdata/manuel/sync/wholeCell/20150505/contactDetection/';
pathBP = '/zdata/manuel/sync/wholeCell/20150505skeletonUpdate/bpc/';
pathGC = '/zdata/manuel/sync/wholeCell/20150505skeletonUpdate/gcl/';
filesBP = dir([pathBP '*.nml']);
filesGC = dir([pathGC '*.nml']);

for i=1:length(filesBP)
    for j=1:length(filesGC)
        idx = sub2ind([length(filesBP) length(filesGC)], i, j);
        functionH{idx} = @contactDetection;
        inputCell{idx} = {[pathBP filesBP(i).name] [pathGC filesGC(j).name] [outputDir filesBP(i).name(1:end-4) 'TO' filesGC(j).name(1:end-4) '.nml']};
    end
end

startCPU(functionH, inputCell, 'retina contact detection full');

end

