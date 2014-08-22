outputDir = '/zdata/manuel/sync/wholeCell/contactDetection/3rd/';
pathBP = '/zdata/manuel/sync/fromLap/ekSkel/3rdTry/bp/';
pathGC = '/zdata/manuel/sync/fromLap/ekSkel/3rdTry/gcl/';
filesBP = dir([pathBP '*.nml']);
filesGC = dir([pathGC '*.nml']);

display('----------------------------------------------------------------------------------------');
tall = tic;
for i=1:length(filesBP)
    for j=1:length(filesGC)
	tSingle = tic;
	display(['Evaluating contacts from ' filesBP(i).name  ' to ' filesGC(j).name]);
	[skel{i,j} data{i,j}] = contactDetection({[pathBP filesBP(i).name] [pathGC filesGC(j).name]});
	if length(skel{i,j}) > 2
		writeNmlOld([outputDir filesBP(i).name(1:end-4) 'TO' filesGC(j).name(1:end-4) '.nml'], skel{i,j});
	end
	toc(tSingle);
	display('----------------------------------------------------------------------------------------');
	end
end
toc(tall);
%save([outputDir 'data.mat'], 'data', 'skel', '-v7.3');


