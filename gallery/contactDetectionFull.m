function contactDetectionFull()

pathSkel = '/zdata/manuel/sync/fromLap/ekSkel/20150804/all/';
outputDir = '/zdata/manuel/sync/wholeCell/cortex/20150804/contactDetection/';

files = dir([pathSkel '*.nml']);
idx = 1;
for i=1:length(files)
    for j=(i+1):length(files)
        functionH{idx} = @contactDetection;
        inputCell{idx} = {[pathSkel files(i).name] [pathSkel files(j).name] [outputDir files(i).name(1:end-4) 'TO' files(j).name(1:end-4) '.nml']};
        idx = idx + 1;
    end
end
startCPU(functionH, inputCell, 'retina contact detection full');

end

