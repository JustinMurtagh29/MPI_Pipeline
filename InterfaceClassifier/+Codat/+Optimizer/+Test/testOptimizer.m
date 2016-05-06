function tests = testOptimizer
%Optimize f(x) = x^2 as simple test case.
tests = functiontests(localfunctions);
end

function testFunctionGD(testCase)
%Test gd
grad = testCase.TestData.fD;
gd = Codat.Optimizer.gradientDescent(1,0.9);
gd = gd.init(1);
param = 10;
for i = 1:1000
    currGrad = grad(param);
    [param,gd] = gd.optimize(param,currGrad);
end
errorBound = 1e-20;
verifyEqual(testCase,param < errorBound,true);
end

function testFunctionRMSProp(testCase)
%Test rmsProp
grad = testCase.TestData.fD;
rmsOpt = Codat.Optimizer.rmsProp(1,0.9,0.9);
rmsOpt = rmsOpt.init(1);
param = 10;
for i = 1:1000
    currGrad = grad(param);
    [param,rmsOpt] = rmsOpt.optimize(param,currGrad);
end
errorBound = 1e-20;
verifyEqual(testCase,param < errorBound,true);
end

function testFunctionAdam(testCase)
%Test rmsProp
grad = testCase.TestData.fD;
ad = Codat.Optimizer.adam(0.1);
ad = ad.init(1);
param = 10;
for i = 1:1000
    currGrad = grad(param);
    [param,ad] = ad.optimize(param,currGrad);
end
errorBound = 1e-20;
verifyEqual(testCase,param < errorBound,true);
end

function setupOnce(testCase)
testCase.TestData.fD = @(x)x;
end

