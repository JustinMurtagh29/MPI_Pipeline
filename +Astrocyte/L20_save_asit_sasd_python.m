asit_long = load('~/20190729T143522_sasd-pair-table.mat');
asit = asit_long.asiT;
asitt = struct(asit);
%% Strings to double
asitt.data{2} = double(asitt.data{2}); %type
asitt.data{8} = double(asitt.data{8}); %axonclass
asitt.data{9} = double(asitt.data{9}); %targetclass

%% the asit
savedir = '~/GABA/astrocyte/synapses/asit45_sasd/asit/';
for i=1:numel(asitt.Properties.VariableNames)
    var = asitt.data{i};
    save([savedir, asitt.Properties.VariableNames{i}, '.mat'], 'var');
end

%% the sasd table
sasd = asit_long.sasdT;
sasdt = struct(sasd);

savedir = '~/GABA/astrocyte/synapses/asit45_sasd/sasd/';
for i=1:numel(sasdt.Properties.VariableNames)
    var = sasdt.data{i};
    save([savedir, sasdt.Properties.VariableNames{i}, '.mat'], 'var');
end

%%
%{
Keys:
id
type
preAggloId
postAggloId
shid
synIds
borderIds
axonClass
targetClass
pos
area

Type:
  1   PrimarySpine 
  2   SecondarySpine 
  4   Soma 

Axon Class:
  1   Corticocortical 
  2   Thalamocortical 
  4   Other 

Target Class:
  1   Somata 
  2   ProximalDendrite 
  3   SmoothDendrite 
  4   ApicalDendrite 
  5   AxonInitialSegment 
  6   OtherDendrite 

%}



