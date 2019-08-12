%% save syn for python
syn = load('~/GABA/astrocyte/synapses/partiallySplit/syn.mat');
%syn = syn.syn;
table = syn.syn.synapses;
table.type = double(table.type);
tablet = struct(table);

savedir = '~/GABA/astrocyte/synapses/partiallySplit/syn/';
for i = 1:4
    var = tablet.data{i};
    save([savedir, table.Properties.VariableNames{i}, '.mat'], 'var');
end

%% save conn for python

conn = load('~/GABA/astrocyte/synapses/partiallySplit/conn.mat');
axons = conn.conn.axons;
save('~/GABA/astrocyte/synapses/partiallySplit/connAxons.mat', 'axons');
dendrites = conn.conn.dendrites;
save('~/GABA/astrocyte/synapses/partiallySplit/connDends.mat', 'dendrites');

