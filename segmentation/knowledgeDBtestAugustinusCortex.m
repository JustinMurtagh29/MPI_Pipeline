%% Load all necessary data
%clear all; clc;
load /data/knowledgeDBlastTest//augustinusStack2.mat;
seg = {stackBig};
% skel = readKnossosNml('/data/knowledgeDBlastTest/2012-09-28_ex145_07x2_2.nml');
% toDelete = []; % Collect empty skeletons
% for l=1:size(skel,2)
%     if ~isfield(skel{l}, 'nodesNumDataAll')
%         toDelete = [toDelete l];
%     end
% end
% skel(toDelete) = [];
% writeKnossosNml('/p/denseMinicubes/cortex/2012-09-28_ex145_07x2.nml.clean', skel);
% skel = readKnossosNml('/p/denseMinicubes/cortex/2012-09-28_ex145_07x2.nml.clean');
% skel = switchToLocalCoords( skel, [2000 2000 2000], 0 );
% skel = correctSkeletonsToBBox(skel, 300);
% writeKnossosNml('/p/denseMinicubes/cortex/2012-09-28_ex145_07x2.nml.local', skel);
%skel = readKnossosNml('/p/denseMinicubes/cortex/2012-09-28_ex145_07x2.nml.local');
totalPathLength = getPathLength(skel);
eval = evaluateSeg(seg, skel, 1);

eval.general = eval.general(1,1);
eval.split = eval.split(1,1);
eval.merge = eval.merge(1,1);
eval.nodes = eval.nodes;

cubePos = [2001 2001 2001];
seg = seg{1,1};

%% Where to put data
knowledgeDB = '/home/mberning/Desktop/knowledgeDBcortexLast/';

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
seg = uint16(seg);
% save data
save /data/knowledgeDBcortexLast.mat

%% load data
load /data/knowledgeDBcortexLast.mat

%% Create struct for settings.json
settings.name = '2012-09-28_ex145_07x2';
settings.maxCoordinates = [10240 7168 3456];
settings.voxelSize = [11.24 11.24 28];
settings.dataLayer.segmentation.rootFolder = [settings.name '/segmentation/'];
settings.dataLayer.segmentation.filePrefix = [settings.name];
settings.dataLayer.segmentation.class = class(seg);
settings.dataLayer.segmentation.resolutions = [1];
settings.dataLayer.raw.rootFolder = [settings.name '/raw/'];
settings.dataLayer.raw.filePrefix = [settings.name];
settings.dataLayer.raw.class = class(raw);
settings.dataLayer.raw.resolutions = [1 2 4 8];

%% Cubes schreiben
writeKnossosRoi([knowledgeDB settings.dataLayer.segmentation.rootFolder], settings.dataLayer.segmentation.filePrefix, cubePos, seg, class(seg));

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
            missions(task).start.position = skel{currentSkeleton}.nodes(nodeIdx1,1:3) + cubePos;
            missions(task).start.direction = (skel{currentSkeleton}.nodes(nodeIdx2,1:3) - skel{currentSkeleton}.nodes(nodeIdx1,1:3)) ./ norm(skel{currentSkeleton}.nodes(nodeIdx2,1:3) - skel{currentSkeleton}.nodes(nodeIdx1,1:3));
            missions(task).start.id = objIndices{i}(nodeIdx1);
            trueid = objIndices{i}(nodeIdx2);
            % Some geometry
            obj1 = seg == objIndices{i}(nodeIdx1);
            obj2 = seg == objIndices{i}(nodeIdx2);
            temp1 = bwdist(obj1) < 2;
            temp2 = bwdist(obj2) < 2;
            temp = temp1 + temp2;
            coM = regionprops(temp == 2, 'Centroid', 'Area', 'PixelIdxList');
            [~, idx] = max([coM(:).Area]);
            if size(coM,1) > 1
                warning(['More than one overlap region ' num2str(i)]);
            end
            if size(coM,1) == 0
                warning('No overlap');
            end
            missions(task).errorCenter = round(coM(idx).Centroid([2 1 3])) + cubePos;
%             [x, y, z] = getRegion(missions(task).errorCenter);
%             targetVolume = seg(x, y, z);
%             temp = (unique(targetVolume(targetVolume ~= 0)))';
            obj = zeros(size(seg));
            obj(coM(idx).PixelIdxList) = 1;
            obj = bwdist(obj) < 5;
            temp = unique(seg(obj));
            exclude = [0 490 538 546 561];
            for k=1:length(exclude)
                temp(temp == exclude(k)) = [];
            end    
            for k=1:length(temp)
                missions(task).possibleEnds(k).id = temp(k);
                missions(task).possibleEnds(k).probability = sum(sum(sum(targetVolume == temp(k))))/ numel(targetVolume);
            end
        end
    end
end

%% Visualization

for task=1:length(missions)
    figure('Position', [1601 1 1920 1079], 'Renderer', 'OpenGL');
    % Plot starting position & direction
    localPos = missions(task).start.position - cubePos;
    plot3(localPos(2), localPos(1), localPos(3), 'xr', 'LineWidth', 5);
    hold on;
    direct = skel{task}.nodes(nodeIdx2,1:3) - skel{task}.nodes(nodeIdx1,1:3);
    quiver3(localPos(2), localPos(1), localPos(3), direct(2), direct(1), direct(3), 0, 'k', 'LineWidth', 3);
    % Plot error center
    localCenter = missions(task).errorCenter - cubePos;
    plot3(localCenter(2), localCenter(1), localCenter(3), 'ob', 'LineWidth', 20);
    % Plot start object solid
    obj = seg == missions(task).start.id;
    issf = isosurface(obj, .1);
    k = patch(issf);
    set(k, 'FaceColor', 'g', 'EdgeColor', 'none', 'FaceAlpha', .5);
    % Plot other objects transparent
    for i=1:length(missions(task).possibleEnds)        colors = jet(length(missions(task).possibleEnds));
        obj = seg == missions(task).possibleEnds(i).id;
        issf = isosurface(obj, .1);
        l = patch(issf);
        set(l, 'FaceColor', colors(i,:), 'EdgeColor', 'none', 'FaceAlpha', .15);
    end
    view(3);
    daspect([28 28 11]);
    title(['Task: ' num2str(task)]);
    grid on;
    camlight('headlight');
    lighting phong;
    legend({'start position', 'start direction', 'error center', 'start object'});
    %pause;
    saveas(gcf, ['/home/mberning/task' num2str(task, '%.2i') '.png']);
    saveas(gcf, ['/home/mberning/task' num2str(task, '%.2i') '.fig']);
    close all;
end

%% Write metadata to JSON
addpath(genpath('auxiliary/jsonlab'));
savejson('', settings, 'FileName', [knowledgeDB settings.name '/settings.json']);
files = dir([knowledgeDB settings.name '/missions*.json']);
savejson('', missions, 'FileName', [knowledgeDB settings.name '/missions' num2str(size(files,1) + 1, '%.4i') '.json']);

