function tests = testConvFcts
%TESTCONVFCTS Test conv3, conv3fft, conv3fft_2 outputs
fprintf('[%s] - Running testConvFcts:\n',datestr(now));
tests = functiontests(localfunctions);
end

function testValidMode(testCase)
%Test valid mode
fprintf('[%s] - Testing valid mode.\n',datestr(now));
inputChannel = 2;
outputChannel = 3;
A = randn(10,10,10,inputChannel);
b = rand(outputChannel,1);
for filterSize = 3:7;
    w = randn(filterSize,filterSize,filterSize,inputChannel,outputChannel);
    y1 = Codat.NN.conv3fft(A,w,b,[1,1,1],'valid');
    y2 = Codat.NN.conv3fft_2(A,w,b,[1,1,1],'valid');
    y3 = Codat.NN.conv3(A,w,b,'valid');
    diff = max([max(abs(y1(:) - y2(:))),max(abs(y2(:) - y3(:))),max(abs(y1(:) - y3(:))),]);
    verifyEqual(testCase,1e-6 > diff,true);
end
end

function testSamedMode(testCase)
%Test valid mode
fprintf('[%s] - Testing same mode.\n',datestr(now));
inputChannel = 2;
outputChannel = 3;
A = randn(10,10,10,inputChannel);
b = rand(outputChannel,1);
for filterSize = 3:7;
    w = randn(filterSize,filterSize,filterSize,inputChannel,outputChannel);
    y1 = Codat.NN.conv3fft(A,w,b,[1,1,1],'same');
    y2 = Codat.NN.conv3fft_2(A,w,b,[1,1,1],'same');
    y3 = Codat.NN.conv3(A,w,b,'same');
    diff = max([max(abs(y1(:) - y2(:))),max(abs(y2(:) - y3(:))),max(abs(y1(:) - y3(:))),]);
    verifyEqual(testCase,1e-6 > diff,true);
end
end

function testFullMode(testCase)
%Test valid mode
fprintf('[%s] - Testing full mode.\n',datestr(now));
inputChannel = 2;
outputChannel = 3;
A = randn(10,10,10,inputChannel);
b = rand(outputChannel,1);
for filterSize = 3:7;
    w = randn(filterSize,filterSize,filterSize,inputChannel,outputChannel);
    y1 = Codat.NN.conv3fft(A,w,b,[1,1,1],'full');
    y2 = Codat.NN.conv3fft_2(A,w,b,[1,1,1],'full');
    y3 = Codat.NN.conv3(A,w,b,'full');
    diff = max([max(abs(y1(:) - y2(:))),max(abs(y2(:) - y3(:))),max(abs(y1(:) - y3(:))),]);
    verifyEqual(testCase,1e-6 > diff,true);
end
end

function testValidDRegularMode(testCase)
%Test valid mode
fprintf('[%s] - Testing d-regularity with valid mode.\n',datestr(now));
inputChannel = 2;
outputChannel = 3;
d = [2,2,2];
A = randn(10,10,10,inputChannel);
b = rand(outputChannel,1);
for filterSize = 3:5;
    w = randn(filterSize,filterSize,filterSize,inputChannel,outputChannel);
    y1 = Codat.NN.conv3fft(A,w,b,d,'valid');
    y2 = Codat.NN.conv3fft_2(A,w,b,d,'valid');
    w = Codat.CNN.cnn.sparseKernel(w,d);
    y3 = Codat.NN.conv3(A,w,b,'valid');
    diff = max([max(abs(y1(:) - y2(:))),max(abs(y2(:) - y3(:))),max(abs(y1(:) - y3(:))),]);
    verifyEqual(testCase,1e-6 > diff,true);
end
end

function testSameDRegularMode(testCase)
%Test valid mode
fprintf('[%s] - Testing d-regularity with same mode.\n',datestr(now));
inputChannel = 2;
outputChannel = 3;
d = [2,2,2];
A = randn(10,10,10,inputChannel);
b = rand(outputChannel,1);
for filterSize = 3:5;
    w = randn(filterSize,filterSize,filterSize,inputChannel,outputChannel);
    y1 = Codat.NN.conv3fft(A,w,b,d,'same');
    y2 = Codat.NN.conv3fft_2(A,w,b,d,'same');
    w = Codat.CNN.cnn.sparseKernel(w,d);
    y3 = Codat.NN.conv3(A,w,b,'same');
    diff = max([max(abs(y1(:) - y2(:))),max(abs(y2(:) - y3(:))),max(abs(y1(:) - y3(:))),]);
    verifyEqual(testCase,1e-6 > diff,true);
end
end

function testFullDRegularMode(testCase)
%Test valid mode
fprintf('[%s] - Testing d-regularity with full mode.\n',datestr(now));
inputChannel = 2;
outputChannel = 3;
d = [2,2,2];
A = randn(10,10,10,inputChannel);
b = rand(outputChannel,1);
for filterSize = 3:5;
    w = randn(filterSize,filterSize,filterSize,inputChannel,outputChannel);
    y1 = Codat.NN.conv3fft(A,w,b,d,'full');
    y2 = Codat.NN.conv3fft_2(A,w,b,d,'full');
    w = Codat.CNN.cnn.sparseKernel(w,d);
    y3 = Codat.NN.conv3(A,w,b,'full');
    diff = max([max(abs(y1(:) - y2(:))),max(abs(y2(:) - y3(:))),max(abs(y1(:) - y3(:))),]);
    verifyEqual(testCase,1e-6 > diff,true);
end
end