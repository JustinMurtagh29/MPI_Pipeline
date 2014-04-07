outputDir = '/zdata/manuel/sync/wholeCell/contactDetection/';
pathBP = '/zdata/manuel/sync/fromLap/ekSkel/bpcForPaperFinal/';
pathGC = '/zdata/manuel/sync/fromLap/ekSkel/gcForPaperFinal/';
filesBP = dir([pathBP '*.nml']);
filesGC = dir([pathGC '*.nml']);

tic;
for i=1:length(filesBP)
    for j=1:length(filesGC)
	display(['Evaluating contacts from ' filesBP(i).name  ' to ' filesGC(j).name]);
	[skel(i,j) merger(i,j)] = contactDetection({[pathBP filesBP(i).name] [pathGC filesGC(j).name]});
	writeNml([outputDir filesBP(i).name(1:end-4) 'TO' filesGC(j).name(1:end-4) '.nml'], skel);
    end
end
save([outputDir 'data.mat'], 'skel', 'merger');
toc

