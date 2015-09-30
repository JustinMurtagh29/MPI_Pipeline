function contactDetectionFull()

pathSkel = '/gaba/u/mberning/data/retina/retinaN2skeletons/allNice/';
outputDir = '/gaba/u/mberning/data/retina/retinaN2skeletons/results/contactDetection/';

files = dir([pathSkel '*.nml']);
idx = 1;
for i=1:length(files)
    for j=(i+1):length(files)
        inputCell{idx} = {[pathSkel files(i).name] [pathSkel files(j).name] [outputDir files(i).name(1:end-4) 'TO' files(j).name(1:end-4) '.nml']};
        idx = idx + 1;
    end
end
functionH = @contactDetection;
startCPU(functionH, inputCell, 'retina contact detection full');

end

