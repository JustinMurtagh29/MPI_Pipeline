outputDir = '/zdata/manuel/sync/wholeCell/contactDetection/cortex/';
pathPre = '/zdata/manuel/sync/fromLap/07x2skeletons/axonsMHforPaper/';;
pathPost = '/zdata/manuel/sync/fromLap/07x2skeletons/spinyStellateForPaper/';
filesPre = dir([pathPre '*.nml']);
filesPost = dir([pathPost '*.nml']);

for i=1:length(filesPre)
    for j=1:length(filesPost)
        idx = sub2ind([length(filesPre) length(filesPost)], i, j);
        functionH{idx} = @contactDetectionCortex;
        inputCell{idx} = {p [pathPre filesPre(i).name] [pathPost filesPost(j).name] [outputDir filesPre(i).name(1:end-4) 'TO' filesPost(j).name(1:end-4) '.nml']};
    end
end

startCPU(functionH, inputCell, 'cortex gallery full');

