function tests = testKnossosIO()
    tests = functiontests(localfunctions);
end

function testCubeAlignedIO(test)
   cubeSize = [512, 512, 256];
   cubePos = [1, 1, 1];
   
   doIt(test, cubePos, cubeSize);
end

function testCubeNonAlignedIO(test)
    cubeSize = [512, 512, 256];
    cubePos = [123, 234, 345];
    
    doIt(test, cubePos, cubeSize);
end

function doIt(test, pos, cubeSize)
    % build fake data
    data = uint8(255 .* rand(cubeSize));
    
    % build bounding box
    box = repmat(pos(:), [1, 2]);
    box(:, 2) = box(:, 2) + size(data)' - 1;
    
    wkDir = setupDir();
    writeKnossosRoi(wkDir, 'test', pos, data);
    readData = readKnossosRoi(wkDir, 'test', box);
    
    % check for equality
    verifyEqual(test, data(:), readData(:));
end

function dirName = setupDir()
    while true
        dirName = tempname();
        
        % MATLAB could really use do-while
        if ~exist(dirName, 'dir');
            break;
        end
    end
    
    % create dir
    mkdir(dirName);
end