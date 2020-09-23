% get synapses corresponding to a given axon agglomerate
% skel with synapse locations exported as nml
%{
rootFolder = '/tmpscratch/sahilloo/data/H2_3_v2_U1_SubI/pipelineRun_mr2e_wsmrnet/20190702T100708_agglomeration/';
load(fullfile(rootFolder,'agglos.mat')); % load state of agglos of axons and dendrites

outputFolder = fullfile(rootFolder, 'axonsWithSyns');
mkdir(outputFolder)

m = load(fullfile(p.saveFolder,'globalSynapseScores.mat'));
synScores = m.synScores;
edgeIdx = m.edgeIdx;

m = load(fullfile(p.saveFolder,'globalEdges.mat'));
globalEdges = m.edges;

m = load(fullfile(p.saveFolder,'globalBorder.mat'));
borderCom = m.borderCoM;

segmentMeta = load([p.saveFolder 'segmentMeta.mat'], 'point', 'maxSegId');

graph = load([p.saveFolder 'graph.mat']);

synT = -1.67;
synIdx = L4.Synapses.getSynapsePredictions(synScores,synT);
allSynEdgeIdx = edgeIdx(synIdx);
allPreSynId = globalEdges(allSynEdgeIdx,1);
allPostSynId = globalEdges(allSynEdgeIdx,2);

% write skel with synapses
display('Writing skeletons for debugging the process:');
parameters.experiment.name= p.experimentName;
parameters.scale.x = num2str(p.raw.voxelSize(1));
parameters.scale.y = num2str(p.raw.voxelSize(2));
parameters.scale.z = num2str(p.raw.voxelSize(3));
parameters.offset.x = '0';
parameters.offset.y = '0';
parameters.offset.z = '0';

maxSegId = segmentMeta.maxSegId;
%}
% choose which axon agglos to search
agglos = axonsSorted(1:100);
for i=1:numel(agglos)
    curAgglo = agglos{i};
    % find synaptic edge in this agglo
    synFound = ismember(allPreSynId,curAgglo);
    if any(synFound)
    edgeSynIdx = allSynEdgeIdx(synFound);
    preSynId = allPreSynId(synFound);
    postSynId = allPostSynId(synFound);
    
    skel = writeSkelWithSyn(graph.edges, segmentMeta.point, curAgglo, ['axon_' num2str(i,'%03d')], parameters, ...
                        edgeSynIdx, preSynId, postSynId, borderCom);
    skel.write(fullfile(outputFolder,['axonWithSynapses_' num2str(i,'%03d') '.nml']))
    clear skel
    else
    display(sprintf('No synapse found for axon: %d \n',i));
    end
end

function skel = writeSkelWithSyn(edges, com, cc, treeName, parameters,...
                curSynEdgeIdx, preSynId, postSynId, borderCom)
    % Set colors to be used
    colors = distinguishable_colors(1, [0 0 0; 1 1 1]);
    colors(:,4) = 0;
    
    skel = initializeSkeleton(parameters);
    skel.names{1} = treeName;
    skel.colors{1} = colors;
    % Get information to write to skeleton
    theseCoM = com(cc,:);
    idx = all(ismember(edges,cc),2);
    theseEdgesSegId = edges(idx,:);
    % remove duplicate edges
    theseEdgesSegId = unique(theseEdgesSegId, 'rows');
    theseEdgesNodes = changem(double(theseEdgesSegId), 1:size(theseCoM,1), cc);
    clear idx theseEdgesSegId;
    % Write to structure for writeNml
    skel = skel.addTree(treeName, theseCoM, theseEdgesNodes);
    clear theseCoM theseEdgesNodes;
    synapses = table(num2cell(curSynEdgeIdx),num2cell(preSynId),num2cell(postSynId),'VariableNames',{'edgeIdx','presynId','postsynId'});
    skel = L4.Synapses.synapse2Skel(synapses,borderCom,com,skel);
    
end

function skel = initializeSkeleton(parameters)
    skel = skeleton();
    if nargin > 0 && ~isempty(parameters)
    skel.parameters = parameters;
    else
    % Set parameters
    skel.parameters.experiment.name='H2_3_v2_U1_SubI';
    skel.parameters.scale.x = '11.24';
    skel.parameters.scale.y = '11.24';
    skel.parameters.scale.z = '28';
    skel.parameters.offset.x = '0';
    skel.parameters.offset.y = '0';
    skel.parameters.offset.z = '0';
    end
end

