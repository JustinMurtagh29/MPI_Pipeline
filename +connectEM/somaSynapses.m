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
    posx = cellfun(@(x){ceil(double(borders.borderCoM(graph.borderIdx(x(1)),:))./[32,32,16])},somasynapsesEdgeIdx{idx});
%     find out whether they are in the soma region
    goodsynapses{idx} = cellfun(@(x)thiscube5(x(1),x(2),x(3)),posx);
    %document
    connectEM.generateSkeletonFromAgglo([ones(sum(goodsynapses{idx})-1,1), (2:sum(goodsynapses{idx}))'],cell2mat(cellfun(@(x){double(borders.borderCoM(graph.borderIdx(x(1)),:))},somasynapsesEdgeIdx{idx}(goodsynapses{idx}))),{1:sum(goodsynapses{idx})},{num2str(idx)},'tata',sum(goodsynapses{idx}));
end

% added by BS: get synapses for each soma
p = Gaba.getSegParameters('ex145_ROI2017');
[graph, ~, borderMeta] = Seg.IO.loadGraph(p, false);
maxSomaId = max(cellfun(@max, somaAgglos.somaAgglos(:,1)));
postsynSomaLUT = L4.Agglo.buildLUT(synapses.synapses.postsynId, maxSomaId);
somaSynIdx = cellfun(@(x)setdiff(postsynSomaLUT(x), 0), ...
    somaAgglos.somaAgglos(:,1), 'uni', 0);

synScoresT.synT = -1;
synScoresT.otherSynT = 1.5;
toDiscard = cell(numSomas, 1);
for i = 1:numSomas
    [~, toDiscard{i}] = L4.Soma.somaSynapsePostProcessing( ...
        synapses.synapses, synapses.boutons, somaSynIdx{i}, ...
        synapses.isSpineSyn, graph.synScores, synScoresT );
end
somaSynIdx = cellfun(@(x,y)x(~y), somaSynIdx, toDiscard, 'uni', 0);
somaSynEdgeIdx = cellfun(@(x)cell2mat(synapses.synapses.edgeIdx(x)), ...
    somaSynIdx, 'uni', 0);
for idx = 1:numSomas
    posx = arrayfun(@(x){ceil(double(borders.borderCoM(graph.borderIdx(x),:))./[32,32,16])},somaSynEdgeIdx{idx});
    load(['thiscube5_' num2str(idx)],'thiscube5');
    goodsynapses{idx} = cellfun(@(x)thiscube5(x(1),x(2),x(3)),posx);
    skel = L4.Agglo.agglo2Nml(somaSynEdgeIdx{idx}(goodsynapses{idx}), borderMeta.borderCoM);
end