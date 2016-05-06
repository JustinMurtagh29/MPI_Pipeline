function tests = testBatchNormalization
%TESTCONVFCTS Test conv3, conv3fft, conv3fft_2 outputs
fprintf('[%s] - Running testBatchNormalization:\n',datestr(now));
tests = functiontests(localfunctions);
end

function testBatchNormalizationForward(testCase)
X = randn(10,10,10,3);
beta = permute([0;0;0],[4 3 2 1]);
gamma = permute([1;1;1],[4 3 2 1]);
Y = Codat.NN.batchNormalization(X,beta,gamma);
for i = 1:size(Y,4)
    tmp = Y(:,:,:,i);
    verifyEqual(testCase,abs(mean(tmp(:))) < 1e-10,true);
    verifyEqual(testCase,abs(var(tmp(:)) - 1) < 1e-2,true);
end

beta = permute([3;3;3],[4 3 2 1]);
gamma = permute([5;5;5],[4 3 2 1]);
Y = Codat.NN.batchNormalization(X,beta,gamma);
for i = 1:size(Y,4)
    tmp = Y(:,:,:,i);
    verifyEqual(testCase,abs(mean(tmp(:)) - 3) < 1e-10,true);
    verifyEqual(testCase,abs(var(tmp(:)) - 25) < 5e-2,true);
end

% test prediction
X = 2.*randn(10,10,10,3) + 5;
beta = permute([3;3;3],[4 3 2 1]);
gamma = permute([5;5;5],[4 3 2 1]);
muInf = permute([5;5;5],[4 3 2 1]);
sig2Inf = permute([4;4;4],[4 3 2 1]);
Y = Codat.NN.batchNormalization(X,beta,gamma,muInf,sig2Inf,false);
for i = 1:size(Y,4)
    tmp = Y(:,:,:,i);
    verifyEqual(testCase,abs(mean(tmp(:)) - 3) < 6e-1,true);
    verifyEqual(testCase,abs(var(tmp(:)) - 25) < 5,true);
end
end

function testBatchNormalizationBwd(testCase)
X = randn(5,5,5,2);
beta = permute([0;0],[4 3 2 1]);
gamma = permute([1;1],[4 3 2 1]);
[Y,muB,sig2B] = Codat.NN.batchNormalization(X,beta,gamma);
target = randn(size(Y));
[~,DEDY] = squaredError(Y,target);
[gradIn, gradb, gradg] = Codat.NN.batchNormalizationBwd(DEDY,X,muB,sig2B,gamma);
epsilon = 1e-8;
maxDiff = 1e-4;

%test derivative wrt input
numGradIn = zeros(size(gradIn));
for i = 1:numel(X)
    x = X;
    x(i) = x(i) + epsilon;
    Y = Codat.NN.batchNormalization(x,beta,gamma);
    pos = squaredError(Y,target);
    x = X;
    x(i) = x(i) - epsilon;
    Y = Codat.NN.batchNormalization(x,beta,gamma);
    neg = squaredError(Y,target);
    numGradIn(i) = (pos - neg)./(2*epsilon);
end
diff = max(abs(numGradIn(:) - gradIn(:)));
verifyEqual(testCase,maxDiff > diff, true);

%test derivative wrt to parameters
numGradb = zeros(size(gradb));
for i = 1:length(beta)
    b = beta;
    b(i) = b(i) + epsilon;
    Y = Codat.NN.batchNormalization(x,b,gamma);
    pos = squaredError(Y,target);
    b = beta;
    b(i) = b(i) - epsilon;
    Y = Codat.NN.batchNormalization(x,b,gamma);
    neg = squaredError(Y,target);
    numGradb(i) = (pos - neg)./(2*epsilon);
end
diff = max(abs(numGradb(:) - gradb(:)));
verifyEqual(testCase,maxDiff > diff, true);

numGradg = zeros(size(gradg));
for i = 1:length(beta)
    g = gamma;
    g(i) = g(i) + epsilon;
    Y = Codat.NN.batchNormalization(x,beta,g);
    pos = squaredError(Y,target);
    g = gamma;
    g(i) = g(i) - epsilon;
    Y = Codat.NN.batchNormalization(x,beta,g);
    neg = squaredError(Y,target);
    numGradg(i) = (pos - neg)./(2*epsilon);
end
diff = max(abs(numGradg(:) - gradg(:)));
verifyEqual(testCase,maxDiff > diff, true);
end

function [loss, DEDY] = squaredError(y, t)
    loss = 0.5.*sum((y(:) - t(:)).^2);
    DEDY = y - t;
end
