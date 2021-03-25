% plot ASI areas vs fit values output from Stan optimization
clear;

%% Configuration
rootDir = '/tmpscratch/sahilloo/data/Mk1_F6_JS_SubI_v1/pipelineRun_mr2e_wsmrnet/';
connFile = fullfile(rootDir, 'connectome', 'Connectome_20191227T220548-results_20191227T220548-results-auto-spines-v3_SynapseAgglomerates--20191227T220548-results--20191227T220548-results-auto-spines-v3--v1.mat');

info = Util.runInfo();
Util.showRunInfo(info);

asiRunId = '20210325T114245'; % NOTE: Using proofread sasd pairs

info = Util.runInfo();
Util.showRunInfo(info);

%% Load axon-spine interfaces
[curDir, curAsiFile] = fileparts(connFile);
curAsiFile = sprintf('%s__%s_sasdData.mat', curAsiFile, asiRunId);
curAsiFile = fullfile(curDir, curAsiFile);

m = load(curAsiFile);
dataTable = m.dataTable;

% choose which synapses to analyse
sasdType = 'spine-spine'; %
switch sasdType
    case 'spine-spine'
        synapseSearchList1 = {'Spine','Prim','Second'};
        synapseSearchList2 = {'Spine','Prim','Second'};
    case 'shaft-shaft'
        synapseSearchList1 = {'Shaft'};
        synapseSearchList2 = {'Shaft'};
    case 'spine-shaft'
        synapseSearchList1 = {'Spine','Prim','Second'};
        synapseSearchList2 = {'Shaft'};
    otherwise
        error('Pls specify correct sasd synapse types')
end

% SASD pairs
idxPlot = (contains(dataTable.syn1,synapseSearchList1) & contains(dataTable.syn2,synapseSearchList2)) | ...
          (contains(dataTable.syn1,synapseSearchList2) & contains(dataTable.syn2,synapseSearchList1)) ;
x1 = dataTable(idxPlot,:).asi1;
x2 = dataTable(idxPlot,:).asi2;

d = sort([x1,x2],2,'descend'); % sort asi data
x1 = d(:,1); x2 = d(:,2);
clear d idxPlot
assert(numel(x1) == numel(x2))
sprintf('Found %d SASD %s pairs',numel(x1), sasdType)

%% asiT
asiT = struct;
asiT.areas = [x1, x2];
asiT.lareas = log10(asiT.areas);

%% Configuration from Stan
mix = struct('mean', {}, 'std', {});

i = 1;
mix(i).mean = -0.67; mix(i).sd = 0.3; i = i + 1; %#ok
mix(i).mean = -0.38; mix(i).sd = 0.12; i = i + 1; %#ok
mix(i).mean = 0.19; mix(i).sd = 0.89; i = i + 1; %#ok
clear i;

rng(0)
for i=1:numel(mix)
    curPd = makedist('Lognormal','mu',mix(i).mean,'sigma',mix(i).sd);
    mix(i).pd = curPd;
end

%% plot
rng(0)
fitClass = 'Normal';
colors = Util.getRangedColors(0,1,0,1,0,1, numel(mix));
fig = figure;
fig.Color = 'white';
ax = gca;
hold on

% raw SASD ASI data
curData = asiT.lareas(:); % log10 scale
color = 'k';
% make histogram
hgram = histogram(curData,'BinWidth',0.1,...
    'DisplayStyle','stairs',...
    'LineWidth',2,...
    'EdgeColor',color);
h = histfit(curData,[],'normal');
h(1).EdgeColor = [1,1,1];
h(1).FaceColor = [1,1,1];
h(1).Visible = 'off';
h(2).Color = color;
h(2).LineStyle = '-';
% pd of the histfit
pd = fitdist(curData,fitClass);
rawName = {sprintf('Raw data (n=%d)',numel(curData)), ...
           ' ', sprintf('mu=%.2f, sd=%.2f (n=%d)', pd.mu, pd.sigma, numel(curData))};

countN = numel(curData);

% Stan fits of SASD ASI data
for i=1:numel(mix)
    curData = log10(random(mix(i).pd, countN,1));
    color = colors(i,:);
    % make histogram
%     hgram = histogram(curData,'BinWidth',0.1,...
%         'DisplayStyle','stairs',...
%         'LineWidth',2,...
%         'EdgeColor',color);
    h = histfit(curData,[],'normal');
    h(1).EdgeColor = [1,1,1];
    h(1).FaceColor = [1,1,1];
    h(1).Visible = 'off';
    h(2).Color = color;
    h(2).LineStyle = '-';
    fitName{i} = {' ', sprintf('mu=%.2f, sd=%.2f (k=%d)', mix(i).mean, mix(i).sd, i)};
end

curLeg = legend([rawName, reshape([fitName{:}], 1, '')]);
set(curLeg, 'Box', 'Off', 'Location', 'best');

ax.LineWidth = 2;
xlabel('Asi area [log_{10}(um^2)]')
ylabel('Frequency')
set(gca,'FontSize',10)
title(ax, ...
    {info.filename; info.git_repos{1}.hash; sasdType}, ...
    'FontWeight', 'normal', 'FontSize', 10, 'Interpreter','none');
Util.setPlotDefault(ax)
outfile = fullfile(rootDir,'connectome','figures',sprintf('%s-raw-vs-stan-fit.png',sasdType));
export_fig(outfile,'-q101', '-nocrop','-transparent')
close all
