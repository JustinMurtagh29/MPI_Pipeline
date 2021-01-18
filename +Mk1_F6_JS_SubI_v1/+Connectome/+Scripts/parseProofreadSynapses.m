% After extract synaptic locations from connectome and connectomeMeta and generate nmls for manual proofreading
% use proofread nml to extract correct pairs and their post target types
%{
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

conn = load(connFile);
segmentMeta = load([p.saveFolder 'segmentMeta.mat'], 'point', 'maxSegId');

% doing for sasd pairs
idxOut = arrayfun(@(x) numel(x{1})==2, conn.connectome.synIdx);
synIndices = conn.connectome.synIdx(idxOut); % maybe more than one synIdx per contact
synIndicesLookup = transpose(horzcat(synIndices{:}));
contactArea = conn.connectomeMeta.contactArea(idxOut);
contactAreaLookup = transpose(horzcat(contactArea{:}));
%}
%% Specify which synapses file to parse after proofread
%nmlFile = fullfile(rootDir,'connectome','nmls','proofread',[curConnName,sprintf('-proofreadSynapses-%s-finished.nml', curVer)]);

% merged with proofread annotations
nmlFile = fullfile(rootDir,'connectome','nmls','proofread','MkL4-all-SASD-merged-finished.nml');
sprintf('Extracting sasd pairs in %s', nmlFile)

skel = skeleton(nmlFile);
treeNames = skel.names;

%% extract data per sasd pair
countId = 500;
sprintf('Searching pairs till Id %d', countId)
outTable = table;
outTable.idSASD = reshape(1:countId,'',1);
outTable.keep = false(countId,1);
syn1 = {}; % type
syn2 = {}; % type
synIdx1 = [];
synIdx2 = [];
asi1 = [];
asi2 = [];

fun = @(x) str2double(x.id);

for idSASD = 1:countId
    curId = sprintf('sasd-%04d',idSASD);
    curPair = find(contains(treeNames, curId));
%    assert(numel(curPair)==2)
    if isempty(curPair)
        continue;
    end
    curTree1 = treeNames{curPair(1)};
    curTree2 = treeNames{curPair(2)};

    idPre(idSASD) = fun(regexp(curTree1,'pre-(?<id>\d+)-','names'));
    idPost(idSASD) = fun(regexp(curTree1,'post-(?<id>\d+)-','names'));
    curSynIdx1 = fun(regexp(curTree1,'synIdx-(?<id>\d+)','names')); 
    curSynIdx2 = fun(regexp(curTree2,'synIdx-(?<id>\d+)','names'));

    % get contactArea from synIdx, Can be multiple matches so match row
    idxMatch = ismember(synIndicesLookup, [curSynIdx1, curSynIdx2], 'rows');
    flipFlag = false;

    if ~any(idxMatch)
        idxMatch = ismember(synIndicesLookup, [curSynIdx2, curSynIdx1], 'rows');
        flipFlag = true; % flipped syn1 and syn2
    end

    if ~any(idxMatch)
        error(sprintf('Could not find matching indices for the syns for idSASD %d', idSASD))
    end

    if flipFlag
        asiSynIdx1 = contactAreaLookup(idxMatch,2);
        asiSynIdx2 = contactAreaLookup(idxMatch,1);
    else
        asiSynIdx1 = contactAreaLookup(idxMatch,1);
        asiSynIdx2 = contactAreaLookup(idxMatch,2);
    end
    % save
    synIdx1(idSASD) = curSynIdx1;
    synIdx2(idSASD) = curSynIdx2;
    asi1(idSASD) = asiSynIdx1(1); % two pairs were duplicate 
    asi2(idSASD) = asiSynIdx2(1); 

    if contains(curTree1,'True') & ~contains(curTree1,{'merger','ign'})  & contains(curTree2,'True') & ~contains(curTree2,{'merger','ign'}) 
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
outTable.asi1 = reshape(asi1,'',1);
outTable.asi2 = reshape(asi2,'',1);

% keep only correct pairs
dataTable = outTable(outTable.keep,:); % keep true only

curLimX = [-4, +1];
curLimY = [-4, +1];

% 1. plot spine-spine
fig = figure;
fig.Color = 'white';
ax = gca;
hold on
idxPlot = contains(dataTable.syn1,{'Spine','Prim','Second'}) & contains(dataTable.syn2,{'Spine','Prim','Second'});
x1 = dataTable(idxPlot,:).asi1;
x2 = dataTable(idxPlot,:).asi2;
%scatter(x1,x2,'kx');
scatter(log10(x1),log(x2),'kx');
count = sum(idxPlot);

%[r1, p1] = corr(x1, x2,'rows','complete');
curFit = fitlm(log10(x1), log10(x2));
plot(curLimX(:), curFit.predict(curLimX(:)), 'k--');
curLeg = legend(sprintf('y = %.2f + %.2fx (R² = %.2f)', ...
    curFit.Coefficients.Estimate, curFit.Rsquared.Ordinary));
set(curLeg, 'Box', 'Off', 'Location', 'South');

ax.YAxis.Limits = curLimX;
ax.XAxis.Limits = curLimY;
axis('square')
ax.LineWidth = 2;
xlabel('Asi 1 area [log(um^2)]')
ylabel('Asi 2 area [log(um^2)]')
set(gca,'FontSize',10)
title({sprintf('Spine-Spine: %d', count)})
outfile = fullfile(rootDir,'connectome','figures','sasd-pairs-ASI-log-spine-spine.png')
export_fig(outfile,'-q101', '-nocrop','-transparent')
close all

% 2. plot shaft-shaft
fig = figure;
fig.Color = 'white';
ax = gca;
hold on;
idxPlot = contains(dataTable.syn1,{'Shaft'}) & contains(dataTable.syn2,{'Shaft'});
x1 = dataTable(idxPlot,:).asi1;
x2 = dataTable(idxPlot,:).asi2;
%scatter(x1,x2,'kx')
scatter(log10(x1),log(x2),'kx');
count = sum(idxPlot);

%[r1, p1] = corr(x1, x2,'rows','complete');
curFit = fitlm(log10(x1), log10(x2));
plot(curLimX(:), curFit.predict(curLimX(:)), 'k--');
curLeg = legend(sprintf('y = %.2f + %.2fx (R² = %.2f)', ...
    curFit.Coefficients.Estimate, curFit.Rsquared.Ordinary));
set(curLeg, 'Box', 'Off', 'Location', 'South');

ax.YAxis.Limits = curLimX;
ax.XAxis.Limits = curLimY;
axis('square')
ax.LineWidth = 2;
xlabel('Asi 1 area [log(um^2)]')
ylabel('Asi 2 area [log(um^2)]')
set(gca,'FontSize',10)
title({sprintf('Shaft-Shaft: %d', count)})
outfile = fullfile(rootDir,'connectome','figures','sasd-pairs-ASI-log-shaft-shaft.png')
export_fig(outfile,'-q101', '-nocrop','-transparent')
close all

% 3. plot spine-shaft
fig = figure;
fig.Color = 'white';
ax = gca;
hold on;
idxTemp1 = contains(dataTable.syn1,{'Spine','Prim','Second'}) & contains(dataTable.syn2,{'Shaft'}); % neck out
idxTemp2 = contains(dataTable.syn2,{'Spine','Prim','Second'}) & contains(dataTable.syn1,{'Shaft'}); % neck out
idxPlot = idxTemp1 | idxTemp2;
x1 = dataTable(idxPlot,:).asi1;
x2 = dataTable(idxPlot,:).asi2;
%scatter(x1,x2,'kx')
scatter(log10(x1),log(x2),'kx');
count = sum(idxPlot);

%[r1, p1] = corr(x1, x2,'rows','complete');
curFit = fitlm(log10(x1), log10(x2));
plot(curLimX(:), curFit.predict(curLimX(:)), 'k--');
curLeg = legend(sprintf('y = %.2f + %.2fx (R² = %.2f)', ...
    curFit.Coefficients.Estimate, curFit.Rsquared.Ordinary));
set(curLeg, 'Box', 'Off', 'Location', 'South');

ax.YAxis.Limits = curLimX;
ax.XAxis.Limits = curLimY;
axis('square')
ax.LineWidth = 2;
xlabel('Asi 1 area [log(um^2)]')
ylabel('Asi 2 area [log(um^2)]')
set(gca,'FontSize',10)
title({sprintf('Spine-Shaft: %d', count)})
outfile = fullfile(rootDir,'connectome','figures','sasd-pairs-ASI-log-spine-shaft.png')
export_fig(outfile,'-q101', '-nocrop','-transparent')
close all

% 3. plot soma-shaft
fig = figure;
fig.Color = 'white';
ax = gca;
hold on;
idxTemp1 = contains(dataTable.syn1,{'Soma'}) & contains(dataTable.syn2,{'Shaft'});
idxTemp2 = contains(dataTable.syn2,{'Soma'}) & contains(dataTable.syn1,{'Shaft'});
idxPlot = idxTemp1 | idxTemp2;
x1 = dataTable(idxPlot,:).asi1;
x2 = dataTable(idxPlot,:).asi2;
%scatter(x1,x2,'kx')
scatter(log10(x1),log(x2),'kx');
count = sum(idxPlot);

%[r1, p1] = corr(x1, x2,'rows','complete');
curFit = fitlm(log10(x1), log10(x2));
plot(curLimX(:), curFit.predict(curLimX(:)), 'k--');
curLeg = legend(sprintf('y = %.2f + %.2fx (R² = %.2f)', ...
    curFit.Coefficients.Estimate, curFit.Rsquared.Ordinary));
set(curLeg, 'Box', 'Off', 'Location', 'South');

ax.YAxis.Limits = curLimX;
ax.XAxis.Limits = curLimY;
axis('square')
ax.LineWidth = 2;
xlabel('Asi 1 area [log(um^2)]')
ylabel('Asi 2 area [log(um^2)]')
set(gca,'FontSize',10)
title({sprintf('Soma-Shaft: %d', count)})
outfile = fullfile(rootDir,'connectome','figures','sasd-pairs-ASI-log-soma-shaft.png')
export_fig(outfile,'-q101', '-nocrop','-transparent')
close all


%% statistics
sprintf('False: %d', sum(~outTable.keep))
sprintf('Spine-Spine: %d', sum(contains(dataTable.syn1,'Spine') & contains(dataTable.syn2,'Spine')))
sprintf('Spine-Prim: %d', sum( (contains(dataTable.syn1,'Spine') & contains(dataTable.syn2,'Prim')) | ...
                                (contains(dataTable.syn1,'Prim') & contains(dataTable.syn2,'Spine')) ))
sprintf('Spine-Second: %d', sum( (contains(dataTable.syn1,'Spine') & contains(dataTable.syn2,'Second')) | ...
                                (contains(dataTable.syn1,'Second') & contains(dataTable.syn2,'Spine')) ))
sprintf('Prim-Prim: %d', sum( (contains(dataTable.syn1,'Prim') & contains(dataTable.syn2,'Prim')) | ...
                                (contains(dataTable.syn1,'Prim') & contains(dataTable.syn2,'Prim')) ))

sprintf('Second-Second: %d', sum( (contains(dataTable.syn1,'Second') & contains(dataTable.syn2,'Second')) | ...
                                (contains(dataTable.syn1,'Second') & contains(dataTable.syn2,'Second')) ))

sprintf('Prim-Second: %d', sum( (contains(dataTable.syn1,'Prim') & contains(dataTable.syn2,'Second')) | ...
                                (contains(dataTable.syn1,'Second') & contains(dataTable.syn2,'Prim')) ))

sprintf('Spine-Shaft: %d', sum( (contains(dataTable.syn1,'Spine') & contains(dataTable.syn2,'Shaft')) | ...
                                (contains(dataTable.syn1,'Shaft') & contains(dataTable.syn2,'Spine')) ))
sprintf('Prim-Shaft: %d', sum( (contains(dataTable.syn1,'Prim') & contains(dataTable.syn2,'Shaft')) | ...
                                (contains(dataTable.syn1,'Shaft') & contains(dataTable.syn2,'Prim')) ))
sprintf('Second-Shaft: %d', sum( (contains(dataTable.syn1,'Second') & contains(dataTable.syn2,'Shaft')) | ...
                                (contains(dataTable.syn1,'Shaft') & contains(dataTable.syn2,'Second')) ))
sprintf('Shaft-Shaft: %d', sum( (contains(dataTable.syn1,'Shaft') & contains(dataTable.syn2,'Shaft')) | ...
                                (contains(dataTable.syn1,'Shaft') & contains(dataTable.syn2,'Shaft')) ))

sprintf('Soma-Soma: %d', sum( (contains(dataTable.syn1,'Soma') & contains(dataTable.syn2,'Soma')) | ...
                                (contains(dataTable.syn1,'Soma') & contains(dataTable.syn2,'Soma')) ))

sprintf('Soma-Shaft: %d', sum( (contains(dataTable.syn1,'Soma') & contains(dataTable.syn2,'Shaft')) | ...
                                (contains(dataTable.syn1,'Shaft') & contains(dataTable.syn2,'Soma')) ))

sprintf('Neck-Shaft: %d', sum( (contains(dataTable.syn1,'Neck') & contains(dataTable.syn2,'Shaft')) | ...
                                (contains(dataTable.syn1,'Shaft') & contains(dataTable.syn2,'Neck')) ))

function type = funType(x)
    if isempty(x)
        type = 'Spine';
    else
        type = x.type;
    end
end
