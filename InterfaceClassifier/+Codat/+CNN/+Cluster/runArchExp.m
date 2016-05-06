% run a set of experiments for a specific architecture on the cluster
% specify a "saveFolder" where everything is stored
% optinally specify a "name" under which files are stored

if ~exist('saveFolder','var') || isempty(saveFolder)
    error('Specify saveFolder.');
end
if ~exist('name','var') || isempty(name)
    name = 'exp';
end

% % define architecture
% cLayer = 9;
% featureMaps = [1 24 24 24 24 32 32 32 32 1];
% filterSize = {[7 7 3],[5 5 3],[4 4 2],[5 5 3],[5 5 3],[5 5 2],[5 5 2],[3 3 2],[3 3 2]};
% dropout = 0.25.*ones(10,1);
% stride = {[1 1 1]};
% shortcut = [0 0 0 2 0 4 0 6 0 0];
% batchNorm = false;
% maxPool = [false, false, false, true, false, false, false, false, false];
% nonLinearity = 'tanh';
% loss = 'squared';

% cLayer = 14;
% featureMaps = [1 16 16 16 24 24 24 24 32 32 32 32 32 32 1];
% filterSize = {[7 7 3],[5 5 3],[5 5 3],[4 4 2],[5 5 2],[3 3 2],[3 3 2],[3 3 2],[3 3 1],[3 3 2],[3 3 2],[3 3 2],[3 3 2],[3 3 1]};
% dropout = 0.25.*ones(15,1);
% stride = {[1 1 1]};
% batchNorm = false;
% maxPool = false(15,1);
% maxPool(10) = true;
% shortcut = [0 0 0 2 0 4 0 6 0 8 0 10 0 12 0];
% nonLinearity = 'tanh';
% loss = 'squared';

cLayer = 23;
featureMaps = [1 16 16 16 16 16 16 16 24 24 24 24 24 24 24 32 32 32 32 32 32 32 32 1];
filterSize = {[7 7 2],[3 3 2],[3 3 1],[3 3 2],[3 3 2],[3 3 2],[3 3 2],[3 3 2],[3 3 1],[3 3 2],[3 3 2],[3 3 2],[3 3 2],[3 3 2],[3 3 2],[3 3 2],[3 3 2],[3 3 2],[3 3 2],[3 3 2],[3 3 2],[3 3 2],[3 3 1]};
shortcut = [0 0 0 2 0 4 0 6 0 8 0 10 0 12 0 14 0 16 0 18 0 20 0 0];
batchNorm = false;
dropout = 0.25.*ones(24,1);
stride = {[1 1 1]};
maxPool = false;
nonLinearity = 'tanh';
loss = 'squared';

optimizer = Codat.Optimizer.rmsProp(1e-5,0.9);

%starting cnet (only weights are important)
cnet = Codat.CNN.cnn(cLayer, featureMaps, filterSize, stride, ...
                         maxPool, dropout, shortcut, batchNorm, ...
                         nonLinearity, loss, optimizer);
W = cnet.W;
cnets = cell(0);

% %plain architecture
% cnets{end + 1} = Codat.CNN.cnn(cLayer, featureMaps, filterSize, stride, ...
%                          maxPool, 0.*dropout, 0.*shortcut, false, ...
%                          nonLinearity, loss, optimizer);
% cnets{end + 1} = Codat.CNN.cnn(cLayer, featureMaps, filterSize, stride, ...
%                          maxPool, 0.*dropout, 0.*shortcut, false, ...
%                          nonLinearity, loss, optimizer);
% cnets{end - 1}.W = W;
% cnets{end}.W = W;
% cnets{end - 1}.optimizer.learningRate = 1e-4;
% cnets{end}.optimizer.learningRate = 1e-5;

%dropout arch
cnets{end + 1} = Codat.CNN.cnn(cLayer, featureMaps, filterSize, stride, ...
                         maxPool, dropout, 0.*shortcut, false, ...
                         nonLinearity, loss, optimizer);
cnets{end + 1} = Codat.CNN.cnn(cLayer, featureMaps, filterSize, stride, ...
                         maxPool, dropout, 0.*shortcut, false, ...
                         nonLinearity, loss, optimizer);
cnets{end - 1}.W = W;
cnets{end}.W = W;
cnets{end - 1}.optimizer.learningRate = 1e-4;
cnets{end}.optimizer.learningRate = 1e-5;

% %bn arch
% cnets{end + 1} = Codat.CNN.cnn(cLayer, featureMaps, filterSize, stride, ...
%                          maxPool, 0.*dropout, 0.*shortcut, batchNorm, ...
%                          nonLinearity, loss, optimizer);
% cnets{end + 1} = Codat.CNN.cnn(cLayer, featureMaps, filterSize, stride, ...
%                          maxPool, 0.*dropout, 0.*shortcut, batchNorm, ...
%                          nonLinearity, loss, optimizer);
% cnets{end - 1}.W = W;
% cnets{end}.W = W;
% cnets{end - 1}.optimizer.learningRate = 1e-2;
% cnets{end}.optimizer.learningRate = 1e-3;

% %bn + dropout arch
% cnets{end + 1} = Codat.CNN.cnn(cLayer, featureMaps, filterSize, stride, ...
%                          maxPool, dropout, 0.*shortcut, batchNorm, ...
%                          nonLinearity, loss, optimizer);
% cnets{end + 1} = Codat.CNN.cnn(cLayer, featureMaps, filterSize, stride, ...
%                          maxPool, dropout, 0.*shortcut, batchNorm, ...
%                          nonLinearity, loss, optimizer);
% cnets{end - 1}.W = W;
% cnets{end}.W = W;
% cnets{end - 1}.optimizer.learningRate = 1e-2;
% cnets{end}.optimizer.learningRate = 1e-3;

% % residual learning arch
% cnets{end + 1} = Codat.CNN.cnn(cLayer, featureMaps, filterSize, stride, ...
%                          maxPool, 0.*dropout, shortcut, false, ...
%                          nonLinearity, loss, optimizer);
% cnets{end + 1} = Codat.CNN.cnn(cLayer, featureMaps, filterSize, stride, ...
%                          maxPool, 0.*dropout, shortcut, false, ...
%                          nonLinearity, loss, optimizer);
% cnets{end - 1}.W = W;
% cnets{end}.W = W;
% cnets{end - 1}.optimizer.learningRate = 1e-3;
% cnets{end}.optimizer.learningRate = 1e-4;

% residual learning arch + drop
cnets{end + 1} = Codat.CNN.cnn(cLayer, featureMaps, filterSize, stride, ...
                         maxPool, dropout, shortcut, false, ...
                         nonLinearity, loss, optimizer);
cnets{end + 1} = Codat.CNN.cnn(cLayer, featureMaps, filterSize, stride, ...
                         maxPool, dropout, shortcut, false, ...
                         nonLinearity, loss, optimizer);
cnets{end - 1}.W = W;
cnets{end}.W = W;
cnets{end - 1}.optimizer.learningRate = 1e-3;
cnets{end}.optimizer.learningRate = 1e-4;

% % residual learning + bn arch
% cnets{end + 1} = Codat.CNN.cnn(cLayer, featureMaps, filterSize, stride, ...
%                          maxPool, 0.*dropout, shortcut, batchNorm, ...
%                          nonLinearity, loss, optimizer);
% cnets{end + 1} = Codat.CNN.cnn(cLayer, featureMaps, filterSize, stride, ...
%                          maxPool, 0.*dropout, shortcut, batchNorm, ...
%                          nonLinearity, loss, optimizer);
% cnets{end - 1}.W = W;
% cnets{end}.W = W;
% cnets{end - 1}.optimizer.learningRate = 1e-2;
% cnets{end}.optimizer.learningRate = 1e-3;

%save reference cnet
save([saveFolder filesep name '_ref.mat'],'cnet');

%run experiments on workers
cluster = getCluster('gpu');
names = { ...
%         [name '_plain_1'], [name '_plain_2'], ...
        [name '_drop_1'], [name '_drop_2'], ...
%          [name '_bn_1'], [name '_bn_2'], ...
%          [name '_dbn_1'], [name '_dbn_2'], ...
         [name '_res_1'], [name '_res_2'], ...
%          [name '_drop_res_1'], [name '_drop_res_2']};
        };
% names = {[name '_plain_1'], [name '_plain_2'], ...
%          [name '_res_1'], [name '_res_2']};
jobs = Codat.CNN.Cluster.train('membraneWeighted',cnets,saveFolder,names,cluster);
clear names cnets cLayer featureMaps filterSize stride maxPool dropout shortcut batchNorm nonLinearity loss optimizer
