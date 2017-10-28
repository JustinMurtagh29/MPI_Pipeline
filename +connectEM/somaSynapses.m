thiscube = readKnossosRoi('/gaba/wKcubes/Connectomics department/2012-09-28_ex145_07x2_ROI2017/segmentation/8/','2012-09-28_ex145_07x2_ROI2016_corrected_mag8', [1,128*12;1,128*12;1,128*8], 'uint32');
somaAgglos = load('/gaba/u/mberning/results/pipeline/20170217_ROI/aggloState/center_somas.mat');
for idx = [1 : 15, 17:27]
    idx
    thiscube2 = ismember(thiscube, somaAgglos.result{idx,2});
    %connectEM.makesurf(thiscube2,[num2str(idx) '.issf']);
    thiscube3 = nlfilter3(thiscube2, @max, [4, 4, 2]);
    thiscube4 = imopen(thiscube3, ones([7,7,7]));
    thiscube5 = imdilate(thiscube4, ones([3,3,3]));
    %connectEM.makesurf(thiscube5,['a' num2str(idx) '.issf']);
    save(['thiscube5_' num2str(idx)],'thiscube5');
end
m = load(['/gaba/u/mberning/results/pipeline/20170217_ROI/allParameter.mat']);
m.p.connectome.synFile = fullfile(p.connectome.saveFolder, 'SynapseAgglos_v2.mat');
synapses = load('/gaba/u/mberning/results/pipeline/20170217_ROI/connectomeState/SynapseAgglos_v2.mat');
borders = load('/gaba/u/mberning/results/pipeline/20170217_ROI/globalBorder.mat');
graph = load('/gaba/u/mberning/results/pipeline/20170217_ROI/graphNewNew.mat','borderIdx');
lookup = repelem(1:height(synapses.synapses),cellfun(@length,synapses.synapses{:,3}));
for idx = [1:15, 17:27]
    idx
    somasynapsesEdgeIdx = synapses.synapses{lookup(ismember(cell2mat(synapses.synapses{:,3}),somaAgglos.result{idx,2})),1};
    load(['thiscube5_' num2str(idx)],'thiscube5');
    posx = cellfun(@(x){ceil(double(borders.borderCoM(graph.borderIdx(x(1)),:))./[32,32,16])},somasynapsesEdgeIdx);
    goodsynapses{idx} = cellfun(@(x)thiscube5(x(1),x(2),x(3)),posx);
end
