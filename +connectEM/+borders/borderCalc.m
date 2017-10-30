function borderCalc(idx)
load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat');
load('/gaba/u/mberning/results/pipeline/20170217_ROI/aggloState/dendrites_09.mat');
load('/gaba/u/mberning/results/pipeline/20170217_ROI/aggloState/axons_08_a.mat');
seg = readKnossosRoi('/gaba/u/mberning/results/pipeline/20170217_ROI/globalSeg/','2012-09-28_ex145_07x2_ROI2016_corrected_mag1', p.local(idx).bboxSmall, 'uint32');
[edges, ind] = connectEM.borders.codeBenedikt(seg);


dendriteNodes = cat(1, dendrites(indBigDends).nodes);
dendriteIDs = dendriteNodes(:,4);

axonNodes = cat(1, axons(indBigAxons).nodes);
axonIDs = axonNodes(:,4);

[~, findAxon] = ismember(edges, axonIDs);
[~, findDendrite] = ismember(edges, dendriteIDs);
interClassEdges = (findAxon(:,1)>0 & findDendrite(:, 2)>0) | (findAxon(:,2)>0 & findDendrite(:, 1)>0);

findingsA = [edges(interClassEdges,:) ind(interClassEdges)];
bigAxons=axons(indBigAxons);
findAxon2 = findAxon;
axonIDs(end+1) = NaN;
findAxon2(findAxon2==0) = length(axonIDs);
findAxon2 = axonIDs(findAxon2);
for idx2 = 1 : length(bigAxons)
    if mod(idx2,100)==0
        idx2
    end
    findingsB{idx2} = ind(any(ismember(findAxon2,bigAxons(idx2).nodes(:,4)),2));
end
save(['/tmpscratch/kboerg/borders/borders_' num2str(idx)],'findingsA','findingsB');

