%% load connectome
conn = load('\tmpscratch\sahilloo\data\Mk1_F6_JS_SubI_v1\pipelineRun_mr2e_wsmrnet\connectome\Connectome_20191227T220548-results_20191227T220548-results-auto-spines-v3_SynapseAgglomerates--20191227T220548-results--20191227T220548-results-auto-spines-v3--v1.mat')

%% plot sparse connectome
i = conn.connectome.edges(:,1);
j = conn.connectome.edges(:,2);
s = cellfun(@numel, conn.connectome.synIdx);
m = size(conn.axons,1); % row dim, pre
n = size(conn.dendrites,1); % col dim, post
connMat = sparse(i,j,s,m,n);
spy(connMat)
colorbar
xlabel('Postsyn dendrites')
ylabel('Presyn axons')
title('Sparse contactome from connectome data')

%%
maxSyn = 12; % exclude outliers
Util.log('Excluding syn > %d', maxSyn)
ss = s;
ss(ss>maxSyn) = maxSyn;
data = sparse(i,j,ss,m,n);

fig = figure();
fig.Color = 'white';
ax = axes(fig);

colorMap = colormap(jet);
colorMap(1,:) = [0,0,0];
climits = [0.0001, maxSyn];

[ii, jj, Mnnz] = find(data); %// get nonzero values and its positions
scatter(1,1,0.5,0) %// make sure the first color corresponds to 0 value.
hold on
scatter(ii,jj,0.5,Mnnz); %// do the actual plot of the nonzero values

% daspect([1 1 1])
colormap(colorMap)
set(gca,'CLim',climits, 'TickDir','out')
colorbar
xlabel('Postsyn dendrites')
ylabel('Presyn axons')
title('Sparse connectome')

saveas(gcf,fullfile(conn.info.param.rootDir,'connectome','sparseConnectome.png'))

%% plot sparse contactome
i = conn.contactome.aggloIdx(:,1);
j = conn.contactome.aggloIdx(:,2);
s = cellfun(@numel,conn.contactome.edgeIdx);
m = size(conn.axons,1); % row dim, pre
n = size(conn.dendrites,1); % col dim, post
connMat = sparse(i,j,s,m,n);

maxSyn = max(s); % exclude outliers
Util.log('Excluding syn > %d', maxSyn)
ss = s;
ss(ss>maxSyn) = maxSyn;
data = sparse(i,j,ss,m,n);

fig = figure();
fig.Color = 'white';
ax = axes(fig);

colorMap = colormap(jet);
%colorMap(1,:) = [0,0,0];
%climits = [0.0001, maxSyn];

[ii, jj, Mnnz] = find(data); %// get nonzero values and its positions
scatter(1,1,0.5,0) %// make sure the first color corresponds to 0 value.
hold on
scatter(ii,jj,0.5,Mnnz); %// do the actual plot of the nonzero values

% daspect([1 1 1])
colormap(colorMap)
%set(gca,'CLim',climits, 'TickDir','out')
colorbar
xlabel('Postsyn dendrites')
ylabel('Presyn axons')
title('Sparse contactome')

saveas(gcf,fullfile(conn.info.param.rootDir,'connectome','sparseContactome.png'))


