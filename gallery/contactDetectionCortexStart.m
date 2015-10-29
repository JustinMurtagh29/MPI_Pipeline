function contactDetectionCortexStart(p)

outputDir = ['/gaba/u/mberning/results/contactDetection/' datestr(clock,30) '/'];
pathPre = '/gaba/u/mberning/data/cortex/07x2skel/axonsForDC/';;
pathPost = '/gaba/u/mberning/data/cortex/07x2skel/ssForDC/';
filesPre = dir([pathPre '*.nml']);
filesPost = dir([pathPost '*.nml']);

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

for i=1:length(filesPre)
    for j=1:length(filesPost)
        idx = sub2ind([length(filesPre) length(filesPost)], i, j);
        inputCell{idx} = {p [pathPre filesPre(i).name] [pathPost filesPost(j).name] [outputDir filesPre(i).name(1:end-4) 'TO' filesPost(j).name(1:end-4) '.nml']};
    end
end

functionH = @contactDetectionCortex;
startCPU(functionH, inputCell, 'cortex gallery full');

end

