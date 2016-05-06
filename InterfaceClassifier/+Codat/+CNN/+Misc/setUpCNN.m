function cnet = setUpCNN( arch, opt, noOutputChannels)
%SETUPCNN Return some predefined CNN architectures.
% INPUT arch: Architecture name (see below).
%       opt: Optimizer name (Standard sgd).
%       loss: Loss function.
%       noOutputChannels: Integer specifying the number of output channels
%             (Options. Default is 1).
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if ~exist('opt','var') || isempty(opt)
    opt = 'sgd';
    disp('Setting optimizer to sgd.');
end
if ~exist('noOutputChannels','var')
    noOutputChannels = [];
end

switch arch
    case 'test'
        cLayer = 3;
        featureMaps = [1 10 10 1];
        filterSize = {[3 3 3],[3 3 3],[1 1 1]};
        dropout = [0.1 0.25 0.25 0];
        stride = {[2 2 2],[1 1 1],[1 1 1]};
        maxPool = [false, false, false];
        batchNorm = false;
        nonLinearity = {'tanh','tanh','tanh'};
        loss = 'squared';

    case 'shallow'
        cLayer = 3;
        featureMaps = [1 16 16 1];
        filterSize = {[8 8 4],[7 7 3],[7 7 3]};
        dropout = [0.2 0.2 0.2 0];
        stride = {[1 1 1]};
        maxPool = [true,true,false];
        batchNorm = false;
        nonLinearity = 'tanh';
        loss = 'squared';

    case '2d'
        cLayer = 9;
        featureMaps = [1 64 64 64 64 64 64 128 128 1];
        filterSize = {[5 5 1],[5 5 1],[5 5 1],[5 5 1],[5 5 1],[5 5 1],[3 3 1],[3 3 1],[3 3 1]};
        dropout = [0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2];
        stride = {[1 1 1],[1 1 1],[2 2 2],[1 1 1],[1 1 1],[2 2 2],[1 1 1],[1 1 1],[1 1 1]};
        maxPool = false;
        batchNorm = false;
        nonLinearity = 'tanh';
        loss = 'squared';

    case 'cnn1'
        cLayer = 9;
        featureMaps = [1 16 16 16 16 16 16 32 32 1];
        filterSize = {[6 6 4],[5 5 2],[5 5 2],[5 5 2],[5 5 2],[5 5 2],[3 3 2],[3 3 2],[3 3 2]};
        dropout = [0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2];
        stride = {[1 1 1],[1 1 1],[1 1 1],[1 1 1],[1 1 1],[1 1 1],[1 1 1],[1 1 1],[1 1 1]};
        maxPool = [false,false,true,false,false,true,false,false,false];
        batchNorm = false;
        nonLinearity = 'tanh';
        loss = 'squared';

    case 'cnn1_rce'
        cLayer = 9;
        featureMaps = [1 16 16 16 16 16 16 32 32 1];
        filterSize = {[6 6 4],[5 5 2],[5 5 2],[5 5 2],[5 5 2],[5 5 2],[3 3 2],[3 3 2],[3 3 2]};
        dropout = [0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2];
        stride = {[1 1 1],[1 1 1],[1 1 1],[1 1 1],[1 1 1],[1 1 1],[1 1 1],[1 1 1],[1 1 1]};
        maxPool = [false,false,true,false,false,true,false,false,false];
        batchNorm = false;
        nonLinearity = {'relu','relu','relu','relu','relu','relu','relu','relu','sigmoid'};
        loss = 'cross-entropy';

    case 'cnn1_large'
        cLayer = 9;
        featureMaps = [1 24 24 24 32 32 32 32 64 1];
        filterSize = {[6 6 4],[5 5 2],[5 5 2],[5 5 2],[5 5 2],[5 5 2],[3 3 2],[3 3 2],[3 3 2]};
        dropout = [0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2];
        stride = {[1 1 1],[1 1 1],[1 1 1],[1 1 1],[1 1 1],[1 1 1],[1 1 1],[1 1 1],[1 1 1]};
        maxPool = [false,false,true,false,false,true,false,false,false];
        batchNorm = false;
        nonLinearity = 'tanh';
        loss = 'squared';

    case 'cnn1_larger_rce'
        cLayer = 9;
        featureMaps = [1 24 24 24 32 32 32 32 64 1];
        filterSize = {[6 6 4],[5 5 2],[5 5 2],[5 5 2],[5 5 2],[5 5 2],[3 3 2],[3 3 2],[3 3 2]};
        dropout = [0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2];
        stride = {[1 1 1],[1 1 1],[1 1 1],[1 1 1],[1 1 1],[1 1 1],[1 1 1],[1 1 1],[1 1 1]};
        maxPool = [false,false,true,false,false,true,false,false,false];
        batchNorm = false;
        nonLinearity = {'relu','relu','relu','relu','relu','relu','relu','relu','sigmoid'};
        loss = 'cross-entropy';

    case 'cnn2'
        cLayer = 9;
        featureMaps = [1 16 16 24 24 24 32 32 32 1];
        filterSize = {[7 7 3],[5 5 3],[4 4 2],[5 5 3],[5 5 3],[5 5 2],[5 5 2],[3 3 2],[3 3 2]};
        dropout = 0.25.*ones(10,1);
        stride = {[1 1 1]};
        shortcut = [0 0 0 0 0 0 0 0 0 0];
        batchNorm = false;
        maxPool = [false, false, false, true, false, false, false, false, false];
        nonLinearity = 'tanh';
        loss = 'squared';

    case 'cnn2_large'
        cLayer = 9;
        featureMaps = [1 24 24 24 32 32 32 64 64 1];
        filterSize = {[7 7 3],[5 5 3],[4 4 2],[5 5 3],[5 5 3],[5 5 2],[5 5 2],[3 3 2],[3 3 2]};
        dropout = 0.25.*ones(10,1);
        stride = {[1 1 1]};
        shortcut = [0 0 0 0 0 0 0 0 0 0];
        batchNorm = false;
        maxPool = [false, false, false, true, false, false, false, false, false];
        nonLinearity = 'tanh';
        loss = 'squared';
        
    case 'cnn2_huge'
        cLayer = 9;
        featureMaps = [1 32 32 32 48 48 48 64 64 1];
        filterSize = {[7 7 3],[5 5 3],[4 4 2],[5 5 3],[5 5 3],[5 5 2],[5 5 2],[3 3 2],[3 3 2]};
        dropout = 0.25.*ones(10,1);
        stride = {[1 1 1]};
        shortcut = [0 0 0 2 0 4 0 6 0 0];
        batchNorm = false;
        maxPool = [false, false, false, true, false, false, false, false, false];
        nonLinearity = 'tanh';
        loss = 'squared';

    case 'cnn2_ece'
        cLayer = 9;
        featureMaps = [1 16 16 24 24 24 32 32 32 1];
        filterSize = {[7 7 3],[5 5 3],[4 4 2],[5 5 3],[5 5 3],[5 5 2],[5 5 2],[3 3 2],[3 3 2]};
        dropout = 0.25.*ones(10,1);
        stride = {[1 1 1]};
        shortcut = [0 0 0 2 0 4 0 6 0 0];
        batchNorm = false;
        maxPool = [false, false, false, true, false, false, false, false, false];
        nonLinearity = {'elu','elu','elu','elu','elu','elu','elu','elu','sigmoid'};
        loss = 'cross-entropy';

    case 'cnn2_ece_large'
        cLayer = 9;
        featureMaps = [1 16 16 16 32 32 32 64 64 1];
        filterSize = {[7 7 3],[5 5 3],[4 4 2],[5 5 3],[5 5 3],[5 5 2],[5 5 2],[3 3 2],[3 3 2]};
        dropout = 0.25.*ones(10,1);
        stride = {[1 1 1]};
        shortcut = [0 0 0 2 0 4 0 6 0 0];
        batchNorm = false;
        maxPool = [false, false, false, true, false, false, false, false, false];
        nonLinearity = {'elu','elu','elu','elu','elu','elu','elu','elu','sigmoid'};
        loss = 'cross-entropy';

    case 'cnn3'
        cLayer = 14;
        featureMaps = [1 16 16 16 24 24 24 24 32 32 32 32 32 32 1];
        filterSize = {[7 7 3],[5 5 3],[5 5 3],[4 4 2],[5 5 2],[3 3 2],[3 3 2],[3 3 2],[3 3 1],[3 3 2],[3 3 2],[3 3 2],[3 3 2],[3 3 1]};
        dropout = 0.25.*ones(15,1);
        stride = {[1 1 1]};
        batchNorm = false;
        maxPool = false(15,1);
        maxPool(10) = true;
        shortcut = [0 0 0 2 0 4 0 6 0 8 0 10 0 12 0];
        nonLinearity = 'tanh';
        loss = 'squared';

    case 'cnn4'
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
        
    case 'cnn4_ece'
        cLayer = 21;
        featureMaps = [1 16 16 16 16 16 16 16 24 24 24 24 24 24 32 32 32 32 32 32 32 1];
        filterSize = {[7 7 2],[5 5 2],[5 5 2],[3 3 2],[3 3 2],[3 3 2],[3 3 2],[3 3 2],[3 3 2],[3 3 2],[3 3 2],[3 3 2],[3 3 2],[3 3 2],[3 3 2],[3 3 2],[3 3 2],[3 3 2],[3 3 2],[3 3 2],[3 3 1]};
        shortcut = [0 0 0 2 0 0 4 0 0 7 0 0 10 0 0 13 0 0 16 0 18 0];
        batchNorm = false;
        dropout = 0.*0.2.*ones(22,1);
        stride = {[1 1 1]};
        maxPool = false;
        nonLinearity = repmat({'elu'},22,1);
        nonLinearity{end} = 'sigmoid';
        loss = 'cross-entropy';
        
    case 'cnn5'
        cLayer = 5;
        featureMaps = [1 8 8 16 16 1];
        filterSize = {[8 8 5],[7 7 4],[7 7 3],[7 7 3],[7 7 3]};
        dropout = [0 0 0 0 0 0];
        stride = {[1 1 1],[1 1 1],[1 1 1],[1 1 1],[1 1 1]};
        maxPool = [false, true, false, false, false];
        batchNorm = false;
        nonLinearity = 'tanh';
        loss = 'squared';


    case 'cnn5l'
        cLayer = 5;
        featureMaps = [1 16 16 32 32 1];
        filterSize = {[8 8 5],[7 7 4],[7 7 3],[7 7 3],[7 7 3]};
        dropout = [0 0 0 0 0 0];
        stride = {[1 1 1],[1 1 1],[1 1 1],[1 1 1],[1 1 1]};
        maxPool = [false, true, false, false, false];
        batchNorm = false;
        nonLinearity = 'tanh';
        loss = 'squared';

    case 'ChallengerDeep'
        cLayer = 14;
        featureMaps = [1 16 16 16 16 16 16 16 16 16 16 16 16 16 1];
        filterSize = {[5 5 3],[5 5 3],[5 5 3],[5 5 3],[5 5 3],[5 5 3],[5 5 3],[5 5 3],[5 5 2],[5 5 2],[5 5 2],[5 5 2],[5 5 2],[5 5 2]};
        dropout = 0.2;
        stride = {[1,1,1]};
        maxPool = false;
        batchNorm = false;
        nonLinearity = 'tanh';
        loss = 'squared';

    case 'ChallengerDeep_rce'
        cLayer = 14;
        featureMaps = [1 16 16 16 16 16 16 16 16 16 16 16 16 16 1];
        filterSize = {[5 5 3],[5 5 3],[5 5 3],[5 5 3],[5 5 3],[5 5 3],[5 5 3],[5 5 3],[5 5 2],[5 5 2],[5 5 2],[5 5 2],[5 5 2],[5 5 2]};
        dropout = 0.2;
        stride = {[1,1,1]};
        maxPool = false;
        batchNorm = false;
        nonLinearity = {'relu','relu','relu','relu','relu','relu','relu',...
                        'relu','relu','relu','relu','relu','relu','sigmoid'};
        loss = 'cross-entropy';

    case 'ChallengerDeep_ece'
        cLayer = 14;
        featureMaps = [1 16 16 16 16 16 16 16 16 16 16 16 16 16 1];
        filterSize = {[5 5 3],[5 5 3],[5 5 3],[5 5 3],[5 5 3],[5 5 3],[5 5 3],[5 5 3],[5 5 2],[5 5 2],[5 5 2],[5 5 2],[5 5 2],[5 5 2]};
        dropout = 0.2;
        stride = {[1,1,1]};
        maxPool = false;
        batchNorm = false;
        nonLinearity = {'elu','elu','elu','elu','elu','elu','elu',...
                        'elu','elu','elu','elu','elu','elu','sigmoid'};
        loss = 'cross-entropy';

    case 'SegEM'
        cLayer = 5;
        featureMaps = [1 10 10 10 10 1];
        filterSize = {[11 11 5],[11 11 5],[11 11 5],[11 11 5],[11 11 5]};
        dropout = [0 0 0 0 0 0];
        stride = {[1 1 1],[1 1 1],[1 1 1],[1 1 1],[1 1 1]};
        maxPool = false;
        batchNorm = false;
        nonLinearity = 'tanh';
        loss = 'squared';

    case 'DS'
        cLayer = 30;
        featureMaps = [1 16 16 16 16 16 16 16 16 16 16 16 16 16 16 16 ...
                         16 16 16 16 16 16 16 16 16 16 16 16 16 16 1];
        filterSize = {[3 3 2], [3 3 2], [3 3 2],[3 3 2], [3 3 1], ...
                      [3 3 2], [3 3 2], [3 3 2],[3 3 2], [3 3 1], ...
                      [3 3 2], [3 3 2], [3 3 2],[3 3 2], [3 3 1], ...
                      [3 3 2], [3 3 2], [3 3 2],[3 3 2], [3 3 1], ...
                      [3 3 2], [3 3 2], [3 3 2],[3 3 2], [3 3 1], ...
                      [3 3 2], [3 3 2], [3 3 2],[3 3 2], [3 3 1]};
        dropout = 0.2;
        maxPool = false;
        batchNorm = false;
        stride = {[1, 1, 1]};
        nonLinearity = 'tanh';
        loss = 'squared';

    case 'DS_ece'
        cLayer = 30;
        featureMaps = [1 16 16 16 16 16 16 16 16 16 16 16 16 16 16 16 ...
                         16 16 16 16 16 16 16 16 16 16 16 16 16 16 1];
        filterSize = {[3 3 2], [3 3 2], [3 3 2],[3 3 2], [3 3 1], ...
                      [3 3 2], [3 3 2], [3 3 2],[3 3 2], [3 3 1], ...
                      [3 3 2], [3 3 2], [3 3 2],[3 3 2], [3 3 1], ...
                      [3 3 2], [3 3 2], [3 3 2],[3 3 2], [3 3 1], ...
                      [3 3 2], [3 3 2], [3 3 2],[3 3 2], [3 3 1], ...
                      [3 3 2], [3 3 2], [3 3 2],[3 3 2], [3 3 1]};
        dropout = 0.2;
        maxPool = false;
        batchNorm = false;
        stride = {[1, 1, 1]};
        nonLinearity = repmat({'relu'},1,30);
        nonLinearity{end} = 'sigmoid';
        loss = 'cross-entropy';

    case 'VDS'
        cLayer = 40;
        featureMaps = [1 16 16 16 16 16 16 16 16 16 16 16 16 16 16 16 ...
                         16 16 16 16 16 16 16 16 16 16 16 16 16 16 16, ...
                         16 16 16 16 16 16 16 16 16 1];
        filterSize = {[3 3 2], [3 3 2], [3 3 2],[3 3 2], [3 3 1], ...
                      [3 3 2], [3 3 2], [3 3 2],[3 3 2], [3 3 1], ...
                      [3 3 2], [3 3 2], [3 3 2],[3 3 2], [3 3 1], ...
                      [3 3 2], [3 3 2], [3 3 2],[3 3 2], [3 3 1], ...
                      [3 3 2], [3 3 2], [3 3 2],[3 3 2], [3 3 1], ...
                      [3 3 2], [3 3 2], [3 3 2],[3 3 2], [3 3 1], ...
                      [3 3 2], [3 3 2], [3 3 2],[3 3 2], [3 3 1], ...
                      [3 3 2], [3 3 2], [3 3 2],[3 3 2], [3 3 1]};
        dropout = 0.2;
        maxPool = false;
        stride = {[1, 1, 1]};
        nonLinearity = 'tanh';
        loss = 'squared';

    otherwise
        error('Architecture not implemented');
end

if ~isempty(noOutputChannels)
    featureMaps(end) = noOutputChannels;
end

switch opt
    case 'sgd'
        optimizer = Codat.Optimizer.gradientDescent(1e-5,0.9);
    case 'rmsprop'
        optimizer = Codat.Optimizer.rmsProp(1e-5,0.9);
    case 'adam'
        optimizer = Codat.Optimizer.adam(1e-5);
    otherwise
        disp('opt not found. Setting sgd as optimizer');
        optimizer = Codat.Optimizer.gradientDescent(1e-5,0.9);
end

if ~exist('shortcut','var') || isempty(shortcut)
    shortcut = 0;
end
if ~exist('batchNorm','var') || isempty(batchNorm)
    batchNorm = true;
end

cnet = Codat.CNN.cnn(cLayer, featureMaps, filterSize, stride, maxPool, dropout, shortcut, batchNorm, nonLinearity, loss, optimizer);

end
