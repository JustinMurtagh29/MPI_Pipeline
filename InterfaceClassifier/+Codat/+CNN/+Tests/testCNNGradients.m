%Test suite for Codat.CNN
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

function tests = testCNNGradients
%Perform gradient tests for different error and activation funtions.
tests = functiontests(localfunctions);
end

function testTanhActivations(testCase)
% Test squared error function with tanh activations.
disp('Testing squared error with tanh activations.');
cnet = testCase.TestData.cnet;
maxDiff = 5e-4;
diff = cnet.numericalGradientCheck();
verifyEqual(testCase,maxDiff > diff,true);
end

function testReluActivations(testCase)
% Test squared error function with relu activations
disp('Testing squared error with relu activations.');
cnet = testCase.TestData.cnet;
cnet.lossFunction = 'squared';
cnet.nonLinearity = repmat({@Codat.CNN.cnn.relu},cnet.layer,1);
cnet.nonLinearityD = repmat({@Codat.CNN.cnn.reluD},cnet.layer,1);
cnet.b = cellfun(@(x) x + 1,cnet.b,'UniformOutput',false);
diff = cnet.numericalGradientCheck();
maxDiff = 5e-4;
verifyEqual(testCase,maxDiff > diff,true);
end

function testSigmoidActivation(testCase)
% Test squared error function with sigmoid activations.
disp('Testing squared error with sigmoid activations.');
cnet = testCase.TestData.cnet;
cnet.lossFunction = 'squared';
cnet.nonLinearity = repmat({@Codat.CNN.cnn.sigmoid},cnet.layer,1);
cnet.nonLinearityD = repmat({@Codat.CNN.cnn.sigmoidD},cnet.layer,1);
diff = cnet.numericalGradientCheck();
maxDiff = 5e-4;
verifyEqual(testCase,maxDiff > diff,true);
end

function testEluActivations(testCase)
% Test squared error function with elu activations.
disp('Testing squared error with elu activations.');
cnet = testCase.TestData.cnet;
cnet.lossFunction = 'squared';
cnet.nonLinearity = repmat({@Codat.CNN.cnn.elu},cnet.layer,1);
cnet.nonLinearityD = repmat({@Codat.CNN.cnn.eluD},cnet.layer,1);
diff = cnet.numericalGradientCheck();
maxDiff = 5e-4;
verifyEqual(testCase,maxDiff > diff,true);
end

function testCrossEntropyLoss(testCase)
% Test cross-entropy error function with sigmoid activations.
disp('Testing cross-entropy error with sigmoid activations.');
cnet = testCase.TestData.cnet;
cnet.lossFunction = 'cross-entropy';
cnet.nonLinearity = repmat({@Codat.CNN.cnn.sigmoid},cnet.layer,1);
cnet.nonLinearityD = repmat({@Codat.CNN.cnn.sigmoidD},cnet.layer,1);
diff = cnet.numericalGradientCheck();
maxDiff = 5e-4;
verifyEqual(testCase,maxDiff > diff,true);
end

function testSoftmaxLoss(testCase)
% Test softmax error function with softmax activation in last layer.
disp('Testing softmax error with softmax activations.');
cLayer = 3;
featureMaps = [1 2 2 3];
filterSize = {[3 3 3],[3 3 3],[1 1 1]};
stride = {[2 2 2],[1 1 1],[1 1 1]};
dropout = [0.1 0.5 0.5 0.5];
maxPool = [false, true, false];
batchNorm = true;
nonLinearity = {'sigmoid','sigmoid','softmax'};
loss = 'softmax';
optimizer = Codat.Optimizer.gradientDescent(1e-8, 0);
cnet = Codat.CNN.cnn(cLayer, featureMaps, filterSize, stride, maxPool, dropout, 0, batchNorm, nonLinearity, loss, optimizer);
cnet = cnet.setParamsToActvtClass(@double);
diff = cnet.numericalGradientCheck();
maxDiff = 5e-4;
verifyEqual(testCase,maxDiff > diff,true);
end

function testMultipleInputs(testCase)
% Test squared error for multiple input channels.
disp('Testing squared error for multiple input channels.');
cLayer = 3;
featureMaps = [3 2 2 1];
filterSize = {[3 3 3],[3 3 3],[1 1 1]};
stride = {[2 2 2],[1 1 1],[1 1 1]};
maxPool = [false, true, false];
dropout = [0.1 0.5 0.5 0.5];
batchNorm = true;
nonLinearity = 'tanh';
loss = 'squared';
optimizer = Codat.Optimizer.gradientDescent(1e-8, 0);
cnet = Codat.CNN.cnn(cLayer, featureMaps, filterSize, stride, maxPool, dropout, 0, batchNorm, nonLinearity, loss, optimizer);
cnet = cnet.setParamsToActvtClass(@double);
diff = cnet.numericalGradientCheck();
maxDiff = 5e-4;
verifyEqual(testCase,maxDiff > diff,true);
end

function testMultipleOutputs(testCase)
% Test squared error for multiple output channels.
disp('Testing squared error for multiple output channels.');
cLayer = 3;
featureMaps = [1 2 2 3];
filterSize = {[3 3 3],[3 3 3],[1 1 1]};
stride = {[2 2 2],[1 1 1],[1 1 1]};
maxPool = [false, true, false];
dropout = [0.1 0.5 0.5 0.5];
batchNorm = true;
nonLinearity = 'tanh';
loss = 'squared';
optimizer = Codat.Optimizer.gradientDescent(1e-8, 0);
cnet = Codat.CNN.cnn(cLayer, featureMaps, filterSize, stride, maxPool, dropout, 0, batchNorm, nonLinearity, loss, optimizer);
cnet = cnet.setParamsToActvtClass(@double);
diff = cnet.numericalGradientCheck();
maxDiff = 5e-4;
verifyEqual(testCase,maxDiff > diff,true);
end

function testParam2Vec(testCase)
% Test param 2 vec functionality
disp('Test param 2 vec functionality.');
cnet = testCase.TestData.cnet;
cnet.b = cellfun(@(x)randn(size(x)),cnet.b,'UniformOutput',false);
vec = cnet.param2Vec();
cnet2 = cnet.vec2Param(vec);
weightsEqual = all(cellfun(@(x,y)isequal(x,y),cnet.W,cnet2.W));
biasEqual = all(cellfun(@(x,y)isequal(x,y),cnet.b,cnet2.b));
betaEqual = all(cellfun(@(x,y)isequal(x,y),cnet.bn_beta,cnet2.bn_beta));
gammaEqual = all(cellfun(@(x,y)isequal(x,y),cnet.bn_gamma,cnet2.bn_gamma));
verifyEqual(testCase,weightsEqual & biasEqual & betaEqual & gammaEqual,true);
end

function testTranslationInvariance(testCase)
%verify translation invariance
disp('Test translation invariance.');
cnet = testCase.TestData.cnet;
cnet.bn_beta = cellfun(@(x)randn(size(x)),cnet.bn_beta,'UniformOutput',false);
cnet.bn_gamma = cellfun(@(x)randn(size(x)),cnet.bn_gamma,'UniformOutput',false);
cnet.isTraining = false;
A = randn(20,20,20);
A1 = A(1:19,1:19,1:19);
A2 = A(2:20,2:20,2:20);
out1 = cnet.forwardPass(A1);
out2 = cnet.forwardPass(A2);
out1 = out1{end};
out2 = out2{end};
out1 = out1(2:end,2:end,2:end);
out2 = out2(1:end-1,1:end-1,1:end-1);
diff = max(abs(out1(:) - out2(:)));
maxDiff = 1e-10;
verifyEqual(testCase,maxDiff > diff,true);
end

function testFFT2Mode(testCase)
% Test fft1 mode
disp('Test fft2 mode.');
cnet = testCase.TestData.cnet;
cnet = cnet.setConvMode('fft2');
diff = cnet.numericalGradientCheck();
maxDiff = 5e-4;
verifyEqual(testCase,maxDiff > diff,true);
end

function testConvnMode(testCase)
% Test convn mode
disp('Test convn mode.');
cnet = testCase.TestData.cnet;
cnet = cnet.setConvMode('convn');
diff = cnet.numericalGradientCheck();
maxDiff = 5e-4;
verifyEqual(testCase,maxDiff > diff,true);
end

function testModeConsistency(testCase)
% Test convAlgs give equal result
disp('Test convAlgs are equal.');
cnet = testCase.TestData.cnet;
input = randn([5, 5, 5] + cnet.border);
cnet = cnet.setConvMode('fft1');
out1 = cnet.predict(input);
cnet = cnet.setConvMode('fft2');
out2 = cnet.predict(input);
cnet = cnet.setConvMode('convn');
out3 = cnet.predict(input);
diff = max([max(abs(out1(:) - out2(:))),max(abs(out2(:) - out3(:))),max(abs(out1(:) - out3(:)))]);
maxDiff = 1e-14;
verifyEqual(testCase,maxDiff > diff,true);
end

function test2D(testCase)
% Test softmax error function with softmax activation in last layer.
disp('Testing 2d net.');
cLayer = 3;
featureMaps = [1 2 2 3];
filterSize = {[3 3 1],[3 3 1],[1 1 1]};
stride = {[2 2 2],[1 1 1],[1 1 1]};
dropout = [0.1 0.5 0.5 0.5];
maxPool = [false, true, false];
batchNorm = true;
nonLinearity = {'tanh','sigmoid','softmax'};
loss = 'softmax';
optimizer = Codat.Optimizer.gradientDescent(1e-8, 0);
cnet = Codat.CNN.cnn(cLayer, featureMaps, filterSize, stride, maxPool, dropout, 0, batchNorm, nonLinearity, loss, optimizer);
cnet = cnet.setParamsToActvtClass(@double);
diff = cnet.numericalGradientCheck();
maxDiff = 5e-4;
verifyEqual(testCase,maxDiff > diff,true);
end

function testRotTrainingIteration(testCase)
disp('Testing rot training iteration.');
cnet = testCase.TestData.cnet;
maxDiff = 5e-4;
diff = cnet.numericalGradientCheckRot();
verifyEqual(testCase,maxDiff > diff,true);
end

function testShorcutConnections(testCase)
% Test shortcut connections
disp('Testing shortcut connections.');
cLayer = 5;
featureMaps = [1 2 2 2 2 1];
filterSize = {[3 3 3],[3 3 3],[3 3 3],[3 3 3],[3 3 3]};
stride = {[1 1 1]};
dropout = [0.1 0.5 0.5 0.5 0.5 0.5];
shortcut = [0 0 2 0 4 0];
batchNorm = true;
maxPool = false;
nonLinearity = {'tanh','tanh','tanh','tanh','tanh'};
loss = 'squared';
optimizer = Codat.Optimizer.gradientDescent(1e-8, 0);
cnet = Codat.CNN.cnn(cLayer, featureMaps, filterSize, stride, maxPool, dropout, shortcut, batchNorm, nonLinearity, loss, optimizer);
cnet = cnet.setParamsToActvtClass(@double);
diff = cnet.numericalGradientCheck();
maxDiff = 5e-4;
verifyEqual(testCase,maxDiff > diff,true);
end

function testForwardAndPredictEqual(testCase)
% Test if forwardPass and predict give the same output.
disp('Testing shortcut connections.');
cLayer = 5;
featureMaps = [1 2 2 2 2 1];
filterSize = {[3 3 3],[3 3 3],[3 3 3],[3 3 3],[3 3 3]};
stride = {[1 1 1]};
dropout = [0.1 0.5 0.5 0.5 0.5 0.5];
shortcut = [0 0 2 0 4 0];
maxPool = false;
batchNorm = true;
nonLinearity = {'tanh','tanh','tanh','tanh','tanh'};
loss = 'squared';
optimizer = Codat.Optimizer.gradientDescent(1e-8, 0);
cnet = Codat.CNN.cnn(cLayer, featureMaps, filterSize, stride, maxPool, dropout, shortcut, batchNorm, nonLinearity, loss, optimizer);
cnet.isTraining = false;
A = randn(20,20,20,'single');
act = cnet.forwardPass(A);
act = act{end};
pred = cnet.predict(A);
diff = max(abs(pred(:) - act(:)));
maxDiff = 5e-4;
verifyEqual(testCase,maxDiff > diff,true);
end


function setup(testCase)
cLayer = 4;
featureMaps = [1 2 2 2 1];
filterSize = {[3 3 3],[3 3 3],[1,1,1],[1 1 1]};
dropout = [0.1 0.5 0.5 0.5 0.5];
stride = {[2 2 2],[1 1 1],[1 1 1],[1 1 1]};
maxPool = [false, true, false, false];
batchNorm = true;
nonLinearity = 'tanh';
loss = 'squared';
optimizer = Codat.Optimizer.gradientDescent(1e-8, 0);
cnet = Codat.CNN.cnn(cLayer, featureMaps, filterSize, stride, maxPool, dropout, 0, batchNorm, nonLinearity, loss, optimizer);
cnet = cnet.setParamsToActvtClass(@double);
testCase.TestData.cnet = cnet;
end

