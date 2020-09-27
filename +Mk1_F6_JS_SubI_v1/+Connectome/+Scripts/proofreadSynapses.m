% extract synaptic locations from connectome and connectomeMeta and generate nmls for manual proofreading

%% Configuration
rootDir = '/tmpscratch/sahilloo/data/Mk1_F6_JS_SubI_v1/pipelineRun_mr2e_wsmrnet/';
connFile = fullfile(rootDir, 'connectome', 'Connectome_20191227T220548-results_20191227T220548-results-auto-spines-v3_SynapseAgglomerates--20191227T220548-results--20191227T220548-results-auto-spines-v3--v1.mat');

info = Util.runInfo();
Util.showRunInfo(info);

[~, curConnName] = fileparts(connFile);
curVer = 'v1';
experimentName = 'Mk1_F6_JS_SubI_v1_mrnet_wsmrnet_axonDendSpines_ga_20191224T235355optimParams'; % wk dataset with axon_dend_spines agglos as segmentation

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
p = param.p;

conn = load(connFile);
segmentMeta = load([p.saveFolder 'segmentMeta.mat'], 'point', 'maxSegId');

%% Specify which synapses to proofread
sprintf('Extracting sasd pairs:')
idxOut = arrayfun(@(x) numel(x{1})==2, conn.connectome.synIdx);
idPre = conn.connectome.edges(idxOut,1);
idPost = conn.connectome.edges(idxOut,2);
synIdx = conn.connectome.synIdx(idxOut); % cell-array
synComs = conn.connectomeMeta.coms(idxOut); %cell-array

%% Write out nml
outFile = fullfile(rootDir,'connectome','nmls', [curConnName,sprintf('-proofreadSynapses-%s', curVer)]);

sprintf('writing out nmls to %s', outFile)
parameters.experiment.name = experimentName;
parameters.scale.x = num2str(p.raw.voxelSize(1));
parameters.scale.y = num2str(p.raw.voxelSize(2));
parameters.scale.z = num2str(p.raw.voxelSize(3));
parameters.offset.x = '0';
parameters.offset.y = '0';
parameters.offset.z = '0';

maxSegId = segmentMeta.maxSegId;
skel = initializeSkeleton(parameters);

for curSASD = 1:numel(synComs)
    curSynComs = synComs{curSASD};
    curSynIdx = synIdx{curSASD};
    assert(size(curSynComs,1)==numel(curSynIdx))

    for i = 1:size(curSynComs,1) 
        curTreeName = sprintf('sasd-%04d-pre-%d-post-%d-synIdx-%d',curSASD, idPre(curSASD), idPost(curSASD), curSynIdx(i));
        thisCom = curSynComs(i,:);
        skel = skel.addTree(curTreeName, thisCom);
    end
end
skel.write(outFile)

function skel = initializeSkeleton(parameters)
    skel = skeleton();
    if nargin > 0 && ~isempty(parameters)
        skel.parameters = parameters;
    else
    % Set parameters
         skel.parameters.experiment.name='Mk1_F6_JS_SubI_v1_mrnet_wsmrnet';
         skel.parameters.scale.x = '11.24';
         skel.parameters.scale.y = '11.24';
         skel.parameters.scale.z = '28';
         skel.parameters.offset.x = '0';
         skel.parameters.offset.y = '0';
         skel.parameters.offset.z = '0';
    end
end

