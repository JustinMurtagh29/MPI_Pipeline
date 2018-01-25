% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%%
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
outFile = '/home/amotta/Desktop/connectomes.xls';

connFiles = strcat('connectome_axons_', {'16_b'; '17_a'; '18_a'}, '.mat');
connFiles = fullfile(rootDir, 'connectomeState', connFiles);

%% loading data
conns = cellfun(@load, connFiles, 'UniformOutput', false);

%%
vals = cell(0, 2);

for idx = 1:numel(conns)
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
end

vals = transpose(vals);
vals = reshape(vals, 2, [], numel(conns));

out = shiftdim(vals(1, :, 1));
out = cat(2, out, shiftdim(vals(2, :, :)));

t = cell2table(out(2:end, 2:end));
t.Properties.VariableNames = out(1, 2:end);
t.Properties.RowNames = out(2:end, 1);

format('long', 'g'); t