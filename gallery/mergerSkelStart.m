function mergerSkelStart(p)

skelPath = '/gaba/u/mberning/data/cortex/mergerNml/';
outputPath = '/gaba/u/mberning/results/wholeCell/07x2/mergerNml/';

if ~exist(outputPath)
    mkdir(outputPath);
end

files = dir([skelPath '*.nml']);

for i=1:length(files)
    inputCell{i} = {p, skelPath, files(i).name, outputPath};
end

functionH = @mergerSkelTests;
startCPU(functionH, inputCell, 'merger nml test 07x2');

end

