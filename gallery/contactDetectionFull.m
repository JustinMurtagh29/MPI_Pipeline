outputDir = '/zdata/manuel/sync/wholeCell/contactDetection/full2/';
pathBP = '/zdata/manuel/sync/fromLap/ekSkel/full/bpc/';
pathGC = '/zdata/manuel/sync/fromLap/ekSkel/full/gcl/';
filesBP = dir([pathBP '*.nml']);
filesGC = dir([pathGC '*.nml']);

% Weird idea, is something messing up parallel reads
tic;
for i=1:length(filesBP)
    skelBP{i} = parseNml([pathBP filesBP(i).name]);
end
toc;
tic;
for i=1:length(filesGC)
    skelGC{i} = parseNml([pathGC filesGC(i).name]);
end
toc;

for i=1:length(filesBP)
    for j=1:length(filesGC)
        idx = sub2ind([length(filesBP) length(filesGC)], i, j);
        functionH{idx} = @contactDetection;
        inputCell{idx} = {skelBP{i} skelGC{j} [outputDir filesBP(i).name(1:end-4) 'TO' filesGC(j).name(1:end-4) '.nml']};
    end
end

startCPU(functionH(1:1000), inputCell(1:1000), 'cortex galley full');

