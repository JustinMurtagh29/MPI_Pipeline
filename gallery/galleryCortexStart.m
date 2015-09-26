function galleryCortexStart(p)

%skelPath = '/gaba/u/mberning/data/cortex/07x2skel/spinyStellateForDendriteChapter/';
skelPath = '/gaba/u/mberning/data/cortex/07x2skel/misc/';
outputPath = ['/gaba/u/mberning/results/wholeCell/07x2/' datestr(clock,30) '/'];

if ~exist(outputPath)
    mkdir(outputPath);
end

files = dir([skelPath '*.nml']);

for i=1:length(files)
    inputCell{i} = {p, skelPath, files(i).name, outputPath};
end

functionH = @galleryCortex;
startCPU(functionH, inputCell, 'whole cell cortex');

end

