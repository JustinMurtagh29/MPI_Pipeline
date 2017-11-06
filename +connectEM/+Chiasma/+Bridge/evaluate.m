% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
chiasmataFile = fullfile( ...
    rootDir, 'chiasmataSplitting', ...
    '20171104T181213-on-axons-10a-plus-10kE3a', ...
    '20171104T184018_chiasmata.mat');

outputDir = '/home/amotta/Desktop/random-chiasmata';

%% load data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

data = load(chiasmataFile);
axonFile = data.info.param.axonFile;
chiasmaParam = data.info.param.chiasmaParam;
chiasmata = data.chiasmata;
clear data;

axons = load(axonFile);
axonIds = find(axons.indBigAxons);
axons = axons.axons(axonIds);

%% select examples
chiasma = [
    19453, 1;
    158, 2];
chiasma = array2table( ...
    chiasma, 'VariableNames', {'axonId', 'chiasmaId'});

%%
paramForChiasma = transpose(horzcat( ...
    fieldnames(chiasmaParam), struct2cell(chiasmaParam)));
paramForChiasma = Util.modifyStruct(param, paramForChiasma{:});

for curIdx = 1:size(chiasma, 1)
    curAxon = axons(chiasma.axonId(curIdx));
    curChiasmata = chiasmata{chiasma.axonId(curIdx)};
    curChiasmaId = chiasma.chiasmaId(curIdx);
    
    connectEM.Chiasma.Bridge.findBridges( ...
        param, chiasmaParam, curAxon, curChiasmata, curChiasmaId);
end