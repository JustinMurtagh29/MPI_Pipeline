% After extract synaptic locations from connectome and connectomeMeta and generate nmls for manual proofreading
% use proofread nml to extract correct pairs and their post target types

%% Configuration
rootDir = '/tmpscratch/sahilloo/data/Mk1_F6_JS_SubI_v1/pipelineRun_mr2e_wsmrnet/';
connFile = fullfile(rootDir, 'connectome', 'Connectome_20191227T220548-results_20191227T220548-results-auto-spines-v3_SynapseAgglomerates--20191227T220548-results--20191227T220548-results-auto-spines-v3--v1.mat');

info = Util.runInfo();
Util.showRunInfo(info);

[~, curConnName] = fileparts(connFile);
curVer = 'v1';

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
p = param.p;

%% Specify which synapses file to parse after proofread
nmlFile = fullfile(rootDir,'connectome','nmls','proofread',[curConnName,sprintf('-proofreadSynapses-%s-finished.nml', curVer)]);
sprintf('Extracting sasd pairs in %s', nmlFile)

skel = skeleton(nmlFile);
treeNames = skel.names;

%% extract data per sasd pair
countId = 100;
outTable = table;
outTable.idSASD = reshape(1:countId,'',1);
outTable.keep = false(countId,1);
syn1 = {};
syn2 = {};

fun = @(x) str2double(x.id)

for idSASD = 1:countId
    curId = sprintf('sasd-%04d',idSASD);
    curPair = find(contains(treeNames, curId));
    assert(numel(curPair)==2)
    curTree1 = treeNames{curPair(1)};
    curTree2 = treeNames{curPair(2)};

    idPre(idSASD) = fun(regexp(curTree1,'pre-(?<id>\d+)-','names'));
    idPost(idSASD) = fun(regexp(curTree1,'post-(?<id>\d+)-','names'));
    synIdx1(idSASD) = fun(regexp(curTree1,'synIdx-(?<id>\d+)','names')); 
    synIdx2(idSASD) = fun(regexp(curTree2,'synIdx-(?<id>\d+)','names'));
    if contains(curTree1,'True') & contains(curTree2,'True') 
        outTable.keep(idSASD) = true;
        syn1{idSASD} = funType(regexp(curTree1,'(True-(?<type>\w+)','names'));
        syn2{idSASD} = funType(regexp(curTree2,'(True-(?<type>\w+)','names'));
    else
        syn1{idSASD} = 'ign';
        syn2{idSASD} = 'ign';
    end
end

outTable.syn1 = reshape(syn1,'',1);
outTable.syn2 = reshape(syn2,'',1);


%% statistics
sprintf('False: %d', sum(~outTable.keep))

sprintf('Spine-Spine: %d', sum(contains(outTable.syn1,'Spine') & contains(outTable.syn2,'Spine')))
sprintf('Spine-Prim: %d', sum( (contains(outTable.syn1,'Spine') & contains(outTable.syn2,'Prim')) | ...
                                (contains(outTable.syn1,'Prim') & contains(outTable.syn2,'Spine')) ))
sprintf('Spine-Second: %d', sum( (contains(outTable.syn1,'Spine') & contains(outTable.syn2,'Second')) | ...
                                (contains(outTable.syn1,'Second') & contains(outTable.syn2,'Spine')) ))
sprintf('Prim-Prim: %d', sum( (contains(outTable.syn1,'Prim') & contains(outTable.syn2,'Prim')) | ...
                                (contains(outTable.syn1,'Prim') & contains(outTable.syn2,'Prim')) ))

sprintf('Second-Second: %d', sum( (contains(outTable.syn1,'Second') & contains(outTable.syn2,'Second')) | ...
                                (contains(outTable.syn1,'Second') & contains(outTable.syn2,'Second')) ))

sprintf('Prim-Second: %d', sum( (contains(outTable.syn1,'Prim') & contains(outTable.syn2,'Second')) | ...
                                (contains(outTable.syn1,'Second') & contains(outTable.syn2,'Prim')) ))

sprintf('Spine-Shaft: %d', sum( (contains(outTable.syn1,'Spine') & contains(outTable.syn2,'Shaft')) | ...
                                (contains(outTable.syn1,'Shaft') & contains(outTable.syn2,'Spine')) ))
sprintf('Prim-Shaft: %d', sum( (contains(outTable.syn1,'Prim') & contains(outTable.syn2,'Shaft')) | ...
                                (contains(outTable.syn1,'Shaft') & contains(outTable.syn2,'Prim')) ))
sprintf('Second-Shaft: %d', sum( (contains(outTable.syn1,'Second') & contains(outTable.syn2,'Shaft')) | ...
                                (contains(outTable.syn1,'Shaft') & contains(outTable.syn2,'Second')) ))
sprintf('Shaft-Shaft: %d', sum( (contains(outTable.syn1,'Shaft') & contains(outTable.syn2,'Shaft')) | ...
                                (contains(outTable.syn1,'Shaft') & contains(outTable.syn2,'Shaft')) ))

sprintf('Soma-Soma: %d', sum( (contains(outTable.syn1,'Soma') & contains(outTable.syn2,'Soma')) | ...
                                (contains(outTable.syn1,'Soma') & contains(outTable.syn2,'Soma')) ))

sprintf('Soma-Shaft: %d', sum( (contains(outTable.syn1,'Soma') & contains(outTable.syn2,'Shaft')) | ...
                                (contains(outTable.syn1,'Shaft') & contains(outTable.syn2,'Soma')) ))

sprintf('Neck-Shaft: %d', sum( (contains(outTable.syn1,'Neck') & contains(outTable.syn2,'Shaft')) | ...
                                (contains(outTable.syn1,'Shaft') & contains(outTable.syn2,'Neck')) ))







function type = funType(x)
    if isempty(x)
        type = 'Spine';
    else
        type = x.type;
    end
end
