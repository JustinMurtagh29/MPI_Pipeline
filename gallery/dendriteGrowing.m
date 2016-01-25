load('/gaba/u/mberning/skeletons/skelTest.mat');

dendriteIds = gt.segIds(strcmp(gt.label, 'dendrite'));
idx = randi(length(dendriteIds),100,1);
dendriteIds = dendriteIds(idx);
dendriteIds = cellfun(@(x)x(x~=0), dendriteIds, 'UniformOutput', false);

prob = [1:-0.01:.95];
for p_i=1:length(prob)
    for i=1:length(dendriteIds)
        display([num2str(i, '%.5i') '_' num2str(prob(p_i), '%3.2f')]);
        [ids{p_i,i} nrAgglo(p_i,i)] = agglomerateSG_simple(graph, dendriteIds{i}, prob(p_i)); 
    end
end
clear p_i i;

idx = 1;
inputCell = cell(6*100,1);
for p_i=1:length(prob)
    for i=1:length(dendriteIds)
        inputCell{idx} = {p, ids{p_i,i}, ...
            ['/gaba/u/mberning/2012-09-28_ex145_07x2_skeletons/isoTest/dendrite' num2str(i, '%.5i') '_' num2str(prob(p_i), '%3.2f') '.issf' ]};
        idx = idx + 1;
    end
end
functionH = @galleryCortexCCfromSG;
startCPU(functionH, inputCell, 'agglomeration test dendrites');
clear p_i i idx functionH inputCell;

