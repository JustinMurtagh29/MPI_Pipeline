outputDir = '/zdata/manuel/sync/wholeCell/contactDetection/full/';
pathBP = '/zdata/manuel/sync/fromLap/ekSkel/full/bpc/';
pathGC = '/zdata/manuel/sync/fromLap/ekSkel/full/gcl/';
filesBP = dir([pathBP '*.nml']);
filesGC = dir([pathGC '*.nml']);

for i=1:length(filesBP)
    for j=1:length(filesGC)
        idx = sub2ind([length(filesBP) length(filesGC)], i, j);
        functionH{idx} = @contactDetection;
        inputCell{idx} = {[pathBP filesBP(i).name] [pathGC filesGC(j).name] outputDir};
    end
end

startCPU(functionH, inputCell, 'cortex galley full');

