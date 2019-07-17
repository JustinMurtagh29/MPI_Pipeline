asit_long = load('~/GABA/astrocyte/synapses/20190717T144949_synapse-table.mat');
asit = asit_long.synapseTable;
asitt = struct(asit);
%% Strings to double
asitt.data{5} = double(asitt.data{5}); %type
asitt.data{6} = double(asitt.data{6}); %axonclass
asitt.data{7} = double(asitt.data{7}); %targetclass

%%
savedir = '~/GABA/astrocyte/synapses/asit44/';
for i=1:numel(asitt.Properties.VariableNames)
    var = asitt.data{i};
    save([savedir, asitt.Properties.VariableNames{i}, '.mat'], 'var');
end

%%
%{
Keys:
id
area
preAggloId
postAggloId
type
axonClass
targetClass
posPre
posSyn
posPost

Type:
  1   PrimarySpine 
  2   SecondarySpine 
  3   Shaft 
  4   Soma 

Axon Class:
  1   Corticocortical 
  2   Thalamocortical 
  3   Inhibitory 
  4   Other 

Target Class:
  1   ApicalDendrite 
  2   AxonInitialSegment 
  4   OtherDendrite 
  5   SmoothDendrite 
  6   Somata 
  7   WholeCell 

%}



