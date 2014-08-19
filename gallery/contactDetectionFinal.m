outputDir = '/zdata/manuel/sync/wholeCell/contactDetection/_final/';
pathBP = '/zdata/manuel/sync/fromLap/ekSkel/_final/bpc/';
pathGC = '/zdata/manuel/sync/fromLap/ekSkel/_final/gcl/';
filesBP = dir([pathBP '*.nml']);
filesGC = dir([pathGC '*.nml']);

matlabpool 12;
tall = tic;
display('----------------------------------------------------------------------------------------');
parfor i=1:length(filesBP)
    for j=1:length(filesGC)
	display(['Evaluating contacts from ' filesBP(i).name  ' to ' filesGC(j).name]);
	skel = contactDetection({[pathBP filesBP(i).name] [pathGC filesGC(j).name]});
	writeNmlOld([outputDir filesBP(i).name(1:end-4) 'TO' filesGC(j).name(1:end-4) '.nml'], skel);
    end
end
display('----------------------------------------------------------------------------------------');
toc(tall);
matlabpool close;

