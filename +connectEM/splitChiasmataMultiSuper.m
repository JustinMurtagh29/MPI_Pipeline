load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameter.mat', 'p');

temp=load(fullfile(p.saveFolder,'aggloState/axons_06_b.mat'));
temp.axons = temp.axons(temp.indBigAxons);

% get query data, i.e. agglo idx, center idx
queries = connectEM.generateQueriesFromChiasmata('/tmpscratch/kboerg/', temp);

% get taskdata from Heiko
fid = fopen('/tmpscratch/kboerg/IDs_MBKMB_L4_chiasma_axon_queries_26_09_2017.txt');
taskdata = fread(fid);
fclose(fid);
taskdata2 = char(taskdata)';
taskdata3 = strsplit(taskdata2,'\n');
taskdata4 = taskdata3(2:end);
taskdata5 = cellfun(@(x){strsplit(x,',')},taskdata4);
% get FF structure
[ff.segIds, ff.neighbours, ff.filenames, ff.nodes, ff.startNode, ff.comments] = connectEM.lookupNmlMulti(p, {'/tmpscratch/kboerg/L4_chiasma_axon_queries_26_09_2017_nmls/'}, false);
second = @(x)x(2);
ff.filenamesShort = cellfun(@(x)second(strsplit(x,'__')), ff.filenames);

for idx = 1 : length(temp.axons)
    idx
    temp.axons(idx).nodesScaled = bsxfun(@times, temp.axons(idx).nodes(:, 1 : 3), [11.24,11.24,28]);
    % find idx in queries for tasks
    taskids = find(queries(:,1)==idx);
    if isempty(taskids)
        continue
    end
    % find taskids in wK for tasks
    taskids2 = arrayfun(@(x){taskdata5{x}{1}},taskids);
    % find idxs in ff for tasks
    taskids3 = cellfun(@(x)find(strcmp(ff.filenamesShort,x)),taskids2);
    % create task defintions
    tasks = arrayfun(@(x,y)struct('tracings',struct('segids',ff.segIds(x),'nodes',ff.nodes(x)),'centeridx',queries(y,7)),taskids3,taskids);
    % run splitting
    connectEM.splitChiasmataMulti(p,temp.axons(idx),tasks,['/tmpscratch/kboerg/chiasmaSplit26/' num2str(idx) '.mat'])
end