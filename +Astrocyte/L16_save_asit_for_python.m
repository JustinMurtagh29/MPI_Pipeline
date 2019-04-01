asit43_long = load('~/GABA/astrocyte/synapses/asit43_long.mat');
asit = asit43_long.asiT;
asitt = struct(asit);
%% Strings to double
asitt.data{2} = double(asitt.data{2}); %primary spine corresponds to 1
asitt.data{8} = double(asitt.data{8});
asitt.data{9} = double(asitt.data{9});

%%
savedir = '~/GABA/astrocyte/synapses/asit43/';
for i=1:11
    var = asitt.data{i};
    save([savedir, asitt.Properties.VariableNames{i}, '.mat'], 'var');
end

%%
% after the python (notebook 23)
ddict = load('~/GABA/astrocyte/synapses/asit43.mat');

%% save syn for python
syn = load('~/GABA/astrocyte/synapses/syn.mat');
%syn = syn.syn;
table = syn.syn.synapses;
table.type = double(table.type);
tablet = struct(table);

savedir = '~/GABA/astrocyte/synapses/syn/';
for i = 1:4
    var = tablet.data{i};
    save([savedir, table.Properties.VariableNames{i}, '.mat'], 'var');
end