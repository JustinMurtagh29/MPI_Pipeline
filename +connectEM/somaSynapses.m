thiscube = readKnossosRoi(['/gaba/wKcubes/Connectomics department/' ...
    '2012-09-28_ex145_07x2_ROI2017/segmentation/8/'], ...
    '2012-09-28_ex145_07x2_ROI2016_corrected_mag8', ...
    [1,128*12;1,128*12;1,128*8], 'uint32');
somaAgglos = load('/gaba/u/mberning/results/pipeline/20170217_ROI/soma_BS/cells_05.mat');
numSomas = size(somaAgglos.somaAgglos, 1);
for idx = 1 : numSomas
    idx
    thiscube2 = ismember(thiscube, somaAgglos.somaAgglos{idx,1});
    %connectEM.makesurf(thiscube2,[num2str(idx) '.issf']);
    thiscube3 = nlfilter3(thiscube2, @max, [4, 4, 2]);
    thiscube4 = imopen(thiscube3, ones([7,7,7]));
    thiscube5 = imdilate(thiscube4, ones([3,3,3]));
    %connectEM.makesurf(thiscube5,['a' num2str(idx) '.issf']);
    save(['thiscube5_' num2str(idx)],'thiscube5');
end
m = load(['/gaba/u/mberning/results/pipeline/20170217_ROI/allParameter.mat']);
% m.p.connectome.synFile = fullfile(p.connectome.saveFolder, 'SynapseAgglos_v2.mat');
synapses = load('/gaba/u/mberning/results/pipeline/20170217_ROI/connectomeState/SynapseAgglos_v2.mat');
borders = load('/gaba/u/mberning/results/pipeline/20170217_ROI/globalBorder.mat');
graph = load('/gaba/u/mberning/results/pipeline/20170217_ROI/graphNewNew.mat','borderIdx');
lookup = repelem(1:height(synapses.synapses),cellfun(@length,synapses.synapses{:,3}));
somasynapsesEdgeIdx = cell(numSomas, 1);
goodsynapses = cell(numSomas, 1);
for idx = numSomas1:numSomas
    idx
    % get all edges that are part of a synapse
    somasynapsesEdgeIdx{idx} = synapses.synapses{lookup(ismember(cell2mat(synapses.synapses{:,3}),somaAgglos.somaAgglos{idx,1})),1};
    % make sure that there are any
    if isempty(somasynapsesEdgeIdx{idx})
        continue
    end
    % load soma region
    load(['thiscube5_' num2str(idx)],'thiscube5');
    % get the position of all edges as defined in somasynapseEdgeIdx
    posx = cellfun(@(x){ceil(double(borders.borderCoM(graph.borderIdx(x(1)),:))./[32,32,16])},somasynapsesEdgeIdx);
    % find out whether they are in the soma region
    goodsynapses{idx} = cellfun(@(x)thiscube5(x(1),x(2),x(3)),posx);
    %document
    connectEM.generateSkeletonFromAgglo([ones(sum(goodsynapses{idx})-1,1), (2:sum(goodsynapses{idx}))'],cell2mat(cellfun(@(x){double(borders.borderCoM(graph.borderIdx(x(1)),:))},somasynapsesEdgeIdx{idx}(goodsynapses{idx}))),{1:sum(goodsynapses{idx})},{num2str(idx)},'tata',sum(goodsynapses{idx}));
end
