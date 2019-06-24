methodUsed = 'LogitBoost'; %'AdaBoostM1'; % 'LogitBoost';
addpath(genpath('/gaba/u/sahilloo/repos/amotta/matlab/'))

rootDir = '/tmpscratch/sahilloo/data/H2_3_v2_U1_SubI/pipelineRun_mr2e_wsmrnet/';
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;
param.experimentName = 'H2_3_v2_U1_SubI_mr2e_wsmrnet';
import H2_3_v2_SubI.TypeEM.*

info = Util.runInfo();

% load training data
featureSetName = 'segmentAgglomerate';
nmlDir = fullfile(param.saveFolder, ...
     'tracings', 'typeEM');
nmlFiles = fullfile(nmlDir, 'proofread', ...
     {'box-1.nml','box-2.nml', 'box-3.nml', ...
     'box-4.nml','box-5.nml','box-6.nml',...
     'box-7.nml','box-8.nml','box-9.nml'});

% do cross validation
classifiers = cell(0);
results = cell(0);

rng(0);
ids = 1:numel(nmlFiles);
numGlia = [];
numAxon = [];
numDend = [];
numGliaN= [];
numAxonN= [];
numDendN= [];

for curBoxId=1:numel(nmlFiles)
idxTest = curBoxId;
idxTrain = setdiff(ids,curBoxId);

% load test set
    gtTest = TypeEM.GroundTruth.loadSet( ...
            param, featureSetName, nmlFiles(idxTest));
    
    numGlia = vertcat(numGlia, sum(gtTest.label(:,1)>0));
    numAxon = vertcat(numAxon, sum(gtTest.label(:,2)>0));
    numDend = vertcat(numDend, sum(gtTest.label(:,3)>0));
    numGliaN = vertcat(numGliaN, sum(gtTest.label(:,1)<0));
    numAxonN = vertcat(numAxonN, sum(gtTest.label(:,2)<0));
    numDendN = vertcat(numDendN, sum(gtTest.label(:,3)<0));
    
end
Util.save(fullfile(param.saveFolder,'typeEM','crossVal_testBox_statitics.mat'),numGlia, numGliaN, numAxon, numAxonN, numDend, numDendN)
figure
hold on
bar(cat(2,numGlia, numAxon, numDend))
bar(cat(2,numGliaN, numAxonN, numDendN),'FaceColor','none','EdgeColor','g')
legend({'+ glia','+ axon','+ dendrite','- glia','- axon','- dendrite' })
xlabel('Test box Id')
ylabel('Numder of labels')
saveas(gcf,fullfile(param.saveFolder,'typeEM','crossVal_testBox_statistics.png'))
close all

