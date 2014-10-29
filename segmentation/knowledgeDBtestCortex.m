%% Load all necessary data

seg = {stackBig};
skel = readKnossosNml('/p/denseMinicubes/cortex/minicube-januar 2013.042.nml');
toDelete = []; % Collect empty skeletons
for l=1:size(skel,2)
    if ~isfield(skel{l}, 'nodesNumDataAll')
        toDelete = [toDelete l];
    end
end
skel(toDelete) = [];
writeKnossosNml('/p/denseMinicubes/cortex/minicube-januar 2013.042.nml.clean', skel);
skel = readKnossosNml('/p/denseMinicubes/cortex/minicube-januar 2013.042.nml.clean');
skel = switchToLocalCoords( skel, [2000 2000 2000], 0 );
skel = correctSkeletonsToBBox(skel, 300);
writeKnossosNml('/p/denseMinicubes/cortex/minicube-januar 2013.042.nml.local', skel);
skel = readKnossosNml('/p/denseMinicubes/cortex/minicube-januar 2013.042.nml.local');
totalPathLength = getPathLength(skel);
eval = evaluateSeg(seg, skel, 1);

eval.general = eval.general(1,1);
eval.split = eval.split(1,1);
eval.merge = eval.merge(1,1);
eval.nodes = eval.nodes;

cubePos = [2001 2001 2001];
seg = seg{1,1};
% seg = permute(seg{1,1}, [2 1 3]);
% raw = permute(raw, [2 1 3]);

%% Where to put data
knowledgeDB = '/home/mberning/Desktop/knowledgeDBcortex/';

%% Segmentation auswachsen 
sizeCube = [300 300 300];
segTemp = zeros(sizeCube + [2 2 2]);
segTemp(2:sizeCube(1)+1,2:sizeCube(2)+1,2:sizeCube(3)+1) = seg;
for i=2:sizeCube(1)+1
    for j=2:sizeCube(2)+1
        for k=2:sizeCube(3)+1
            if segTemp(i,j,k) == 0
                field = segTemp(i-1:i+1,j-1:j+1,k-1:k+1);
                field = field(:);
                field(field == 0) = [];
                new = mode(field);
                segTemp(i,j,k) = new;
            end
        end
    end
end
segNew = segTemp(2:sizeCube(1)+1,2:sizeCube(2)+1,2:sizeCube(3)+1);
seg = segNew;

% save data
save /data/knowledgeDBcortex2.mat

%% load data
load /data/knowledgeDBcortex2.mat

%% Cubes schreiben
parameter.settings.dataset = '2012-09-28_ex145_07x2';
seg = uint16(seg);

parameter.settings.voxelData.seg.rootFolder = ['voxelData/seg/' parameter.settings.dataset '/mag1/'];
parameter.settings.voxelData.seg.filePrefix = [parameter.settings.dataset '_mag1'];
parameter.settings.voxelData.seg.class = class(seg);
writeKnossosRoi([knowledgeDB parameter.settings.voxelData.seg.rootFolder], parameter.settings.voxelData.seg.filePrefix, cubePos, seg, 'uint16');
parameter.settings.voxelData.raw.rootFolder = ['voxelData/raw/' parameter.settings.dataset '/mag1/'];
parameter.settings.voxelData.raw.filePrefix = [parameter.settings.dataset '_mag1'];
parameter.settings.voxelData.raw.class = class(raw);
writeKnossosRoi([knowledgeDB parameter.settings.voxelData.raw.rootFolder], parameter.settings.voxelData.raw.filePrefix, cubePos, raw);
% parameter.settings.voxelData.flag.rootFolder = ['voxelData/flags/' parameter.settings.dataset '/mag1/'];
% parameter.settings.voxelData.flag.filePrefix = [parameter.settings.dataset '_mag1'];
% parameter.settings.voxelData.flag.class = 'uint8';
% parameter.settings.voxelData.flag.order = {'vesicle' '' '' '' '' '' '' ''};
% writeKnossosRoi([knowledgeDB parameter.settings.voxelData.flag.rootFolder], parameter.settings.voxelData.flag.filePrefix, cubePos, uint8(bw3));

%% Metadata fuer Tasks generieren
task = 0;
for i=1:length(eval.split.idx)
    currentSkeleton = eval.split.idx(i);
    objIndices{i} = seg(sub2ind(size(seg), skel{currentSkeleton}.nodes(:,1), skel{currentSkeleton}.nodes(:,2), skel{currentSkeleton}.nodes(:,3)));
    for j=1:size(skel{eval.split.idx(i)}.edges,1)
        nodeIdx1 = skel{eval.split.idx(i)}.edges(j,1);
        nodeIdx2 = skel{eval.split.idx(i)}.edges(j,2);
        if objIndices{i}(nodeIdx1) ~= objIndices{i}(nodeIdx2) && objIndices{i}(nodeIdx1) ~= 0 && objIndices{i}(nodeIdx2) ~= 0
            task = task + 1;
            parameter.tasks(task).start.position = skel{currentSkeleton}.nodes(nodeIdx1,1:3) + cubePos;
            parameter.tasks(task).start.direction = (skel{currentSkeleton}.nodes(nodeIdx2,1:3) - skel{currentSkeleton}.nodes(nodeIdx1,1:3)) ./ norm(skel{currentSkeleton}.nodes(nodeIdx2,1:3) - skel{currentSkeleton}.nodes(nodeIdx1,1:3));
            parameter.tasks(task).start.id = objIndices{i}(nodeIdx1);
            trueid = objIndices{i}(nodeIdx2);
            % Some geometry
            temp = seg == objIndices{i}(nodeIdx1);
            centerOfMass1 =  regionprops(temp, 'Centroid');
            centerOfMass1 = centerOfMass1.Centroid;
            parameter.tasks(task).start.center =  round(centerOfMass1) + cubePos;
            temp = seg == objIndices{i}(nodeIdx2);
            centerOfMass2 = regionprops(temp, 'Centroid');
            centerOfMass2 = centerOfMass2.Centroid;
%             parameter.tasks(task).end.center = round(centerOfMass2) + cubePos;
            %parameter.tasks(task).errorCenter = round(mean([centerOfMass1; centerOfMass2]));
            centerOfMass2 = round(centerOfMass2);
            [x, y, z] = getRegion(skel{currentSkeleton}.nodes(nodeIdx2,1:3));
            targetVolume = seg(x, y, z);
            temp = (unique(targetVolume(targetVolume ~= 0)))';
            for k=1:length(temp)
                parameter.tasks(task).possibleEnds(k).id = temp(k);
                parameter.tasks(task).possibleEnds(k).probability = sum(sum(sum(targetVolume == temp(k))))/ numel(targetVolume);
            end
        end
    end
end

%% Visualization
% for task=1:length(parameter.tasks)
%     figure('Position', [1601 1 1920 1079], 'Renderer', 'OpenGL');
%     localPos = parameter.tasks(task).start.position - cubePos;
%     plot3(localPos(2), localPos(1), localPos(3), 'xr', 'LineWidth', 5);
%     hold on;
%     quiver3(localPos(2), localPos(1), localPos(3), 10*parameter.tasks(task).start.direction(2), 10*parameter.tasks(task).start.direction(1), 10*parameter.tasks(task).start.direction(3), 5, 'k', 'LineWidth', 2);
%     obj = seg == parameter.tasks(task).start.id;
%     issf = isosurface(obj, .1);
%     k = patch(issf);
%     set(k, 'FaceColor', 'g', 'EdgeColor', 'none', 'FaceAlpha', 1);
%     for i=1:length(parameter.tasks(task).possibleEnds)
%         colors = jet(length(parameter.tasks(task).possibleEnds));
%         obj = seg == parameter.tasks(task).possibleEnds(i).id;
%         issf = isosurface(obj, .1);
%         l = patch(issf);
%         set(l, 'FaceColor', colors(i,:), 'EdgeColor', 'none', 'FaceAlpha', .1);
%     end
%     view(3);
%     daspect([28 28 11]);
%     title(['Task: ' num2str(task)]);
%     grid on;
%     camlight('headlight');
%     lighting phong;
%     %legend({'start position', 'start direction', 'start object', 'end objects', 'Start Object', 'End Object'});
%     for j=1:36
%         camorbit(10,0); camlight('headlight'); pause(.5);
%     end
%     close all;
% end


%% Remove tasks with single possibleEnd
 i=1; 
 while i<=length(parameter.tasks);
     if length(parameter.tasks(i).possibleEnds) == 1
         parameter.tasks(i) = [];
     else
         i = i + 1;
     end
 end 

%% Write metadata to JSON
addpath(genpath('auxiliary/jsonlab'));
savejson('', parameter, [knowledgeDB 'meta.json']);

%% Write series of pngs for Yates
% bw3(bw3 ==1) = 255;
% for i=1:384
%     imwrite(uint8(raw(:,:,i)), [knowledgeDB 'imageSeries/' 'raw' num2str(i, '%.4i') '.png'], 'png');
%     imwrite(uint16(seg(:,:,i)), [knowledgeDB 'imageSeries/' 'seg' num2str(i, '%.4i') '.png'], 'png');
%     imwrite(uint8(bw3(:,:,i)), [knowledgeDB 'imageSeries/' 'flag' num2str(i, '%.4i') '.png'], 'png');
% end

