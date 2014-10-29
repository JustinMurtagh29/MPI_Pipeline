%% Read data from BT_tasks
fileID = fopen('/home/mberning/Desktop/BT_tasks_21-01-2013.csv');
labels = textscan(fileID, [repmat('%s',1,19) '%*[^\n]'], 1, 'Delimiter', ';');
labels = cellfun(@cell2mat, labels, 'UniformOutput', false);
matrix = textscan(fileID, [repmat('%u',1,12) '%s' repmat('%u',1,3) repmat('%s',1,3) '\n'], 'Delimiter', ';');

%% Locate Labels if needed
%whichColumn = find(~cellfun(@isempty, strfind(labels, 'coord')));

%% Look up tasks of tasktype 15
toKeep = matrix{6} == 15 & matrix{12} == 2;
for i=1:length(matrix)
    matrix{i} = matrix{i}(toKeep);
end

%% Find tasks within different buckets belonging to region of interest
bbox =[2001 2100; 2101 2200; 2201 2300];
whichLabel = [9,10,11];
inBbox = cell(3,3);
for i=1:length(whichLabel)
    for j=1:size(bbox,1)
        inBbox{i,j} = bbox(j,1) <= matrix{whichLabel(i)} & bbox(j,2) >= matrix{whichLabel(i)};
    end
end

%% Find tasks and how they fit into ROI
taskArray = zeros(3,3,3);
for x=1:3
    for y=1:3
        for z=1:3
            whichTask = inBbox{1,x} & inBbox{2,y} & inBbox{3,z};
            if any(whichTask)
                taskArray(x,y,z) = matrix{1}(whichTask);
            end
        end
    end
end

%% Stitch together ROI (and renumber segments)
stackDir = '/data/cortexStitch/stackKLEE/';
stackBigWithoutRenumbering = zeros(300,300,300, 'uint32');
stackArray = cell(size(taskArray));
for x=1:3
    for y=1:3
        for z=1:3
            load([stackDir num2str(taskArray(x,y,z)) '.mat']);
            stackArray{x,y,z} = uint32(stack);
            % Where to insert into big stack
            lowerBound = ([x y z] - 1)*100+1;
            upperBound = [x y z]*100;
            % Renumber to remove duplicate IDs
            maxVal = max(max(stackBigWithoutRenumbering(:)),1);
            stackArray{x,y,z} = stackArray{x,y,z} + uint32(stackArray{x,y,z} > 0).*maxVal;
            % Check borders
            stackBigWithoutRenumbering(lowerBound(1):upperBound(1),lowerBound(2):upperBound(2),lowerBound(3):upperBound(3)) = stackArray{x,y,z};
        end
    end
end
stackBig = renumberSegments(stackBigWithoutRenumbering, stackArray, 50);
raw = readKnossosRoi('/z/CortexConnectomics/shared/2012-09-28_ex145_07x2/mag1', '2012-09-28_ex145_07x2_mag1', [2001 2300; 2001 2300; 2001 2300]);