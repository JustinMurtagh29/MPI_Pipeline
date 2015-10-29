function galleryCortexStartFull(p)

load /gaba/u/mberning/data/cortex/07x2skel/all_nodes_solved.mat;
outputPath = ['/gaba/u/mberning/results/wholeCell/07x2/' datestr(clock,30) '/'];

if ~exist(outputPath)
    mkdir(outputPath);
end

idx = 1;
for i=1:length(nodes_export)
    if size(nodes_export{i},1) > 10
        temp{1}.nodes = nodes_export{i}; % to 'mimick' skeleton structure
        inputCell{idx} = {p, temp, [hashes_export{i} '.nml'], outputPath};
        idx = idx + 1;
    end
end

functionH = @galleryCortex;
startCPU(functionH, inputCell, 'whole cell cortex new');

end

