function galleryCortexCCfromSGStart(p)

outputPath = '/zdata/manuel/sync/wholeCell/cortex/CCfromSG/';

if ~exist(outputPath)
    mkdir(outputPath);
end

load([p.saveFolder 'graphNew.mat']);

for i=1:length(graph.ccEdgesJoined.equivalenceClasses)
    functionH{i} = @galleryCortexCCfromSG;;
    inputCell{i} = {p, graph.ccEdgesJoined.equivalenceClasses{i}, [outputPath num2str(i, '%.5i') '.issf']};
end

functionH{1}(inputCell{1}{:});
%startCPU(functionH, inputCell, 'whole cell cortex new');

end

