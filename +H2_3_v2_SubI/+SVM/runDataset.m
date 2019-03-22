% run SVM predictions on the whole dataset
rootDir = '/tmpscratch/sahilloo/data/H2_3_v2_U1_SubI/pipelineRun_mr2e_wsmrnet/';
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

margin = param.tileBorder(:, 2);
assert(isequal(margin, [256; 256; 128]));

info = Util.runInfo();
Util.showRunInfo(info);

% predictions using these trained CNNs
m = load('/gaba/u/bstaffle/data/SyConn/cnn17_enn.mat');
cnn17 = m.cnet;
m = load('/gaba/u/bstaffle/data/SyConn/cnn17_enn.mat');
cnn17_2 = m.cnet;
m = load(['/gaba/u/bstaffle/data/CNN_Training/SVM_CNN/' ...
    'cnn17_SVM_4_large/cnn17_SVM_4_large.mat']);
cnn17_4 = m.cnet;

border = [78, 78, 28]; % cnn17
% prediction options
options.display = 10;
options.gpuDev = false;

%% Run per-cube SVM predictions
taskArgs = arrayfun( ...
    @(local) {local.bboxSmall, local.saveFolder}, ...
    param.local, 'UniformOutput', false);
sharedArgs = {param, border, options, cnn17, cnn17_2, cnn17_4};
sharedArgLocs = [1,4,5,6,7,8];

job = Cluster.startJob( ...
    @H2_3_v2_U1_SubI.SVM.runBox, ...
    taskArgs, 'numOutputs', 0, ...
    'sharedInputs', sharedArgs, ...
    'sharedInputsLocation', sharedArgLocs, ...
    'cluster', {'priority', 100, 'time', '24:00:00', 'memory', 48});

Cluster.waitForJob(job);
Util.log('Finished all box predictions')

