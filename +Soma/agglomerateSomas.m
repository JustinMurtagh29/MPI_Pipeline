load(fullfile('/gaba/u/mberning/results/pipeline/20170217_ROI', 'allParameter.mat'));
graph = load([p.saveFolder 'graphNew.mat'], 'prob', 'edges', 'borderIdx');
meta = load(fullfile(p.saveFolder, 'segmentMeta.mat'), 'segIds', 'point');
load('/gaba/u/mberning/results/pipeline/20170217_ROI/soma/NucleiCoordinates.mat','rp');
gb = load([p.saveFolder 'globalBorder.mat'], 'borderCoM', 'borderSize');

%%choose somaIds
somaIDs = linspace(1,125,125);
%remove glia etc.
somaIDs([6,8,19,21,23,27,28,36,38,41,44,51,57,61,82,87,90,92,95,98,100,104,108,110,117,124]) = [];
somaIDs(61) = []; %makes problems

cells1 = getSomaNodesPar(p, graph, meta, rp, gb, 0.98, 1500, somaIDs);

%second run
somaIDs = [5,29,33,34,39,53,55,65,83,85,97,99,101,113,119];
cells2 = getSomaNodesPar(p, graph, meta, rp, gb, 0.99, 2000, somaIDs);

%third run
somaIDs = [29,33,39,55,83,85,97,99,101,113,119];
cells3 = getSomaNodesPar(p, graph, meta, rp, gb, 0.995, 2250, somaIDs);

save('/tmpscratch/rhesse/soma/results/cells1.mat','cells1')
save('/tmpscratch/rhesse/soma/results/cells2.mat','cells2')
save('/tmpscratch/rhesse/soma/results/cells3.mat','cells3')

%take corrections from run 2&3
somas = cells1;
somas(5,:) = cells2(1,:);
somas(22,:) = [];
somas(25,:) = cells3(2,:);
somas(26,:) = cells2(4,:);
somas(29,:) = [];
somas(39,:) = cells2(6,:);
somas(41,:) = [];
somas(48,:) = cells2(8,:);
somas(64,:) = cells3(5,:);
somas(66,:) = [];
somas(73,:) = cells3(7,:);
somas(74,:) = cells3(8,:);
somas(75,:) = [];
somas(83,:) = [];
somas(87,:) = [];

save('/tmpscratch/rhesse/soma/results/somas_test.mat','somas');
save('/gaba/u/mberning/results/pipeline/20170217_ROI/aggloState/somas_test.mat','somas');

%add somas that merge into another soma
somas(92,:) = cells1(22,:);
somas(93,:) = cells1(43,:);
somas(94,:) = cells1(69,:);
somas(95,:) = cells1(88,:);
somas(96,:) = cells1(93,:);

save('/tmpscratch/rhesse/soma/results/somas_with_merged_somas_test.mat','somas');
save('/gaba/u/mberning/results/pipeline/20170217_ROI/aggloState/somas_with_merged_somas_test.mat','somas');
