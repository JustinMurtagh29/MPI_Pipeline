% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%%
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';

connFiles = { ...
    'connectome_axons-19-a_dendrites-wholeCells-03-v2-classified_spine-syn-clust.mat';
    'connectome_axons-19-a_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat'};
connFiles = fullfile(rootDir, 'connectomeState', connFiles);

info = Util.runInfo();
Util.showRunInfo(info);

%% loading data
conns = cellfun(@load, connFiles, 'UniformOutput', false);

%% collect numbers
connVals = cell(numel(conns), 2);

for idx = 1:numel(conns)
    vals = cell(0, 2);
    conn = conns{idx};
    
   [~, connName] = fileparts(connFiles{idx});
    vals(end + 1, :) = {'Connectome', connName};
    
    vals(end + 1, :) = { ...
        '# axons', numel(conn.axons)};
    vals(end + 1, :) = { ...
        '# postsynaptic targets', numel(conn.dendrites)};
    
    vals(end + 1, :) = { ...
        '# synapses', conn.connectomeMeta.noSynapses};
    vals(end + 1, :) = { ...
        '# spine synapses', sum(conn.axonMeta.spineSynCount)};
    
    vals(end + 1, :) = {...
        '# axons without synapses', ...
        sum(conn.axonMeta.synCount == 0)};
    vals(end + 1, :) = { ...
        '# axons with less than 10 synapses', ...
        sum(conn.axonMeta.synCount < 10)};

    postSynCount = accumarray( ...
        conn.connectome.edges(:, 2), ...
        cellfun(@numel, conn.connectome.synIdx)', ...
       [numel(conn.dendrites), 1]);
    
    vals(end + 1, :) = { ...
        '# targets without synapses', ...
        sum(postSynCount == 0)};
    vals(end + 1, :) = { ...
        '# targets with less than 10 synapses', ...
        sum(postSynCount < 10)};
    
    for classIdx = 1:numel(conn.denClasses)
        vals(end + 1, :) = { ...
            sprintf('# synapses onto %s', conn.denClasses{classIdx}), ...
            sum(conn.classConnectome(:, classIdx))};
    end
    
    connVals{idx, 1} = vals(:, 1);
    connVals{idx, 2} = vals(:, 2);
end

%% show results
rowNames = cat(1, connVals{:, 1});
rowNames = unique(rowNames, 'stable');

tableData = cell(numel(rowNames), numel(conns));
tableData(:) = num2cell(nan);

for idx = 1:numel(conns)
   [~, rows] = ismember(connVals{idx, 1}, rowNames);
    tableData(rows, idx) = connVals{idx, 2};
end

colNames = {'Old', 'New'};
tableData(1, :) = [];
rowNames(1, :) = [];

t = cell2table(tableData);
t.Properties.VariableNames = colNames;
t.Properties.RowNames = rowNames;

format('long', 'g');
disp(t);
