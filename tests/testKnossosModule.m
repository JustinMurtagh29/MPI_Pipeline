function tests = testKnossosModule()
    tests = functiontests(localfunctions);
end

function testMeanAndStd(testCase)
    % prepare cube
    cubeSize = [512, 512, 256];
    cubeElemCount = prod(cubeSize);
    
    myMean = 128;
    myStd = 25;
    
    % build fake raw data
    diff = zeros(cubeSize);
    diff(1:2:cubeElemCount) =  myStd;
    diff(2:2:cubeElemCount) = -myStd;
    
    raw = myMean * ones(cubeSize);
    raw = raw + diff;
    raw = uint8(raw);
    
    % compute mean and std
    [meanVal, stdVal] = ...
        Knossos.calcMeanAndStd(raw);
    
    % use 'eq' to test equality across types
    verifyTrue(testCase, eq(meanVal, myMean));
    verifyTrue(testCase, eq(stdVal, myStd));
end