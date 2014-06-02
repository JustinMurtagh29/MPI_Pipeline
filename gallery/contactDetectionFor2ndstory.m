outputDir = '/zdata/manuel/sync/wholeCell/contactDetection/2nd/';
pathBP = '/zdata/manuel/sync/fromLap/ekSkel/2ndStorySkel/bp/';
pathGC = '/zdata/manuel/sync/fromLap/ekSkel/2ndStorySkel/gcl/';
filesBP = dir([pathBP '*.nml']);
filesGC = dir([pathGC '*.nml']);

tic;
for i=1:length(filesBP)
    for j=1:length(filesGC)
	display(['Evaluating contacts from ' filesBP(i).name  ' to ' filesGC(j).name]);
	[skel{i,j} data{i,j}] = contactDetection({[pathBP filesBP(i).name] [pathGC filesGC(j).name]});
	writeNmlOld([outputDir filesBP(i).name(1:end-4) 'TO' filesGC(j).name(1:end-4) '.nml'], skel{i,j});
    end
end
save([outputDir 'data.mat'], 'data', 'skel', '-v7.3');
toc

