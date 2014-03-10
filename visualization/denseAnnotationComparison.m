%% Take a look at differences between dense annotations
clear all; clc;
directory = 'P:\denseSkeletons\';
files = dir([directory '*.nml']);
for i=1:length(files)
    name{i} = files(i).name(1:end-4);
    skel{i} = parseNml([directory files(i).name]);
    if ~isnumeric(skel{i}{1}.parameters.scale.x)
        skel{i}{1}.parameters.scale.x = str2double(skel{i}{1}.parameters.scale.x);
        skel{i}{1}.parameters.scale.y = str2double(skel{i}{1}.parameters.scale.y);
        skel{i}{1}.parameters.scale.z = str2double(skel{i}{1}.parameters.scale.z);
    end
    scale{i} = skel{i}{1}.parameters.scale;
end
% Fix scale of skeleton 1 Densetracing.nml, somehow 12 by 12 by 24
scale{1} = scale{2};
skel{1}{1}.parameters.scale = scale{1};
% Add bbox of annotated region in voxel
bbox{1} = [1417 1717; 4739 5039; 890 1190];
bbox{2} = [4097 4736; 4481 5248; 2250 2450];
bbox{3} = [3072 3455; 512 895; 2560 2943];
bbox{4} = [1664 2047; 1792 2175; 2176 2559];
for i=1:length(files)
    display(name{i});
    skel{i} = removeEmptySkeletons(skel{i});
    skel{i} = switchToLocalCoords_v2(skel{i}, bbox{i}(:,1)' - [1 1 1]);
    skel{i} = correctSkeletonsToBBox_v2(skel{i}, bbox{i}(:,2)'-bbox{i}(:,1)');
end

%% Calculate volume, path length, number skeletons, number nodes and some aggregate statistics
for i=1:length(name)
    volume(i) = prod((bbox{i}(:,2) - bbox{i}(:,1)) .* [scale{i}.x; scale{i}.y; scale{i}.z])./10^9;
    pathLength(i) = getPathLength(skel{i})./10^3;
    nrSkel(i) = length(skel{i});
    temp = 0;
    for j=1:length(skel{i})
        temp = temp + size(skel{i}{j}.nodes,1);
    end
    nrNodes(i) = temp;
    avgDistNodes(i) = pathLength(i)./nrNodes(i)*10^3;
    pathLengthPerVolume(i) = pathLength(i)./volume(i);
end

%% Plot those statistics
figure;
bar([volume; pathLength; nrSkel; nrNodes; avgDistNodes; pathLengthPerVolume]);
set(gca,'YScale', 'log');
set(gca, 'XTickLabel', {'Volume[microns^3]' 'Path Length[microns]' '# skeletons' '# nodes' 'distance nodes [nm]' 'path length per volume [microns^-2]'});
legend(name);
legend('cortex test', 'cortex training', 'retina test', 'retina training');

%% Get single skeletons out for Moritz
skel{1}{1}.parameters.scale.x = num2str(skel{1}{1}.parameters.scale.x);
skel{1}{1}.parameters.scale.y = num2str(skel{1}{1}.parameters.scale.y);
skel{1}{1}.parameters.scale.z = num2str(skel{1}{1}.parameters.scale.z);
for i=1:length(skel{1})ra
    tempSkel{1}.nodes = skel{1}{i}.nodes;
    tempSkel{1}.nodesAsStruct = skel{1}{i}.nodesAsStruct;
    tempSkel{1}.nodesNumDataAll = skel{1}{i}.nodesNumDataAll;
    tempSkel{1}.parameters = skel{1}{1}.parameters;
    tempSkel{1}.edges = skel{1}{i}.edges;
    tempSkel{1}.thingID = skel{1}{i}.thingID;
    tempSkel{1}.name = skel{1}{i}.name;
    tempSkel{1}.parameters.activeNode.id = tempSkel{1}.nodesAsStruct{1}.id;
    writeNml(['C:\Users\mberning\Desktop\singleSkelDenseAlex\skel' num2str(i, '%.3i') '.nml'], tempSkel);
    clear tempSkel;
end
