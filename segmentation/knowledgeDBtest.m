%% Load all necessary data
map = 5;
algo = 2;
r = 2;
par1 = 4;
par2 = 2;
load([param.dataFolder param.affSubfolder param.affMaps(map).name '.mat']);
load([param.dataFolder param.outputSubfolder param.affMaps(map).name '/seg' num2str(r) '-' num2str(algo) '.mat']);
seg = v{par1,par2};
clear v;
load([param.dataFolder param.outputSubfolder param.affMaps(map).name '/evaluation' num2str(r) '-' num2str(algo) '.mat']);
eval.general = v.general(par1,par2);
eval.split = v.split(par1,par2);
eval.merge = v.merge(par1,par2);
eval.nodes = v.nodes;
clear v;
cubePos = [1665 1793 2177];
seg = permute(seg, [2 1 3]);
raw = permute(raw, [2 1 3]);

%% Where to put data
knowledgeDB = '/home/mberning/Desktop/knowledgeDB/';

%% Prelimiary vesicle detection
rmpath('KLEE');
aff = affX+affY+affZ;
aff = permute(aff, [2 1 3]);

filter1 = -ones(5,5,4);
filter1(2:4,2:4,2:3) = 82 / 18;
% Filter approximately according to histograms
result1 = imfilter(-single(raw), filter1, 'replicate');
bw1 = result1 > 1500;
bw2 = aff > .4;
bw3 = bw2 & bw1;
bw3 = bw3 - bwareaopen(bw3, 40, 6);

% addpath('KLEE');
% KLEE_v4('stack', raw, 'stack_2', bw3);

%% Segmentation auswachsen 
segTemp = zeros(386,386,386);
segTemp(2:385,2:385,2:385) = seg;
for i=2:385
    for j=2:385
        for k=2:385
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
segNew = segTemp(2:385,2:385,2:385);
seg = segNew;

%% save data
save /data/knowledgeDB.mat

%% load data
load /data/knowledgeDB.mat

%% Cubes schreiben
parameter.settings.dataset = 'e_k0563';
seg = uint16(seg);

parameter.settings.voxelData.seg.rootFolder = ['voxelData/seg/' parameter.settings.dataset '/mag1/'];
parameter.settings.voxelData.seg.filePrefix = [parameter.settings.dataset '_mag1'];
parameter.settings.voxelData.seg.class = class(seg);
writeKnossosRoi([knowledgeDB parameter.settings.voxelData.seg.rootFolder], parameter.settings.voxelData.seg.filePrefix, cubePos, seg, 'uint16');
parameter.settings.voxelData.raw.rootFolder = ['voxelData/raw/' parameter.settings.dataset '/mag1/'];
parameter.settings.voxelData.raw.filePrefix = [parameter.settings.dataset '_mag1'];
parameter.settings.voxelData.raw.class = class(raw);
writeKnossosRoi([knowledgeDB parameter.settings.voxelData.raw.rootFolder], parameter.settings.voxelData.raw.filePrefix, cubePos, raw);
parameter.settings.voxelData.flag.rootFolder = ['voxelData/flags/' parameter.settings.dataset '/mag1/'];
parameter.settings.voxelData.flag.filePrefix = [parameter.settings.dataset '_mag1'];
parameter.settings.voxelData.flag.class = 'uint8';
parameter.settings.voxelData.flag.order = {'vesicle' '' '' '' '' '' '' ''};
writeKnossosRoi([knowledgeDB parameter.settings.voxelData.flag.rootFolder], parameter.settings.voxelData.flag.filePrefix, cubePos, uint8(bw3));

%% Metadata fuer Tasks generieren
task = 0;
for i=1:length(eval.split.idx)
    currentSkeleton = eval.split.idx(i);
    objIndices{i} = seg(sub2ind(size(seg), param.skel{currentSkeleton}.nodes(:,1), param.skel{currentSkeleton}.nodes(:,2), param.skel{currentSkeleton}.nodes(:,3)));
    for j=1:size(param.skel{eval.split.idx(i)}.edges,1)
        nodeIdx1 = param.skel{eval.split.idx(i)}.edges(j,1);
        nodeIdx2 = param.skel{eval.split.idx(i)}.edges(j,2);
        if objIndices{i}(nodeIdx1) ~= objIndices{i}(nodeIdx2) && objIndices{i}(nodeIdx1) ~= 0 && objIndices{i}(nodeIdx2) ~= 0
            task = task + 1;
            parameter.tasks(task).start.position = param.skel{currentSkeleton}.nodes(nodeIdx1,1:3) + cubePos;
            parameter.tasks(task).start.direction = (param.skel{currentSkeleton}.nodes(nodeIdx2,1:3) - param.skel{currentSkeleton}.nodes(nodeIdx1,1:3)) ./ norm(param.skel{currentSkeleton}.nodes(nodeIdx2,1:3) - param.skel{currentSkeleton}.nodes(nodeIdx1,1:3));
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
            [x, y, z] = getRegion(centerOfMass2);
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
for task=1:length(parameter.tasks)
    figure('Position', [1601 1 1920 1079]);
    localPos = parameter.tasks(task).start.position - cubePos;
    plot3(parameter.tasks(task).start.center(1), parameter.tasks(task).start.center(2), parameter.tasks(task).start.center(3), 'xr');
    hold on;
    plot3(parameter.tasks(task).end.center(1), parameter.tasks(task).end.center(2), parameter.tasks(task).end.center(3), 'xg');
    plot3(parameter.tasks(task).errorCenter(1), parameter.tasks(task).errorCenter(2), parameter.tasks(task).errorCenter(3), 'xb');
    quiver3(localPos(2), localPos(1), localPos(3), 10*parameter.tasks(task).start.direction(2), 10*parameter.tasks(task).start.direction(1), 10*parameter.tasks(task).start.direction(3), 'y', 'LineWidth', 5);
    obj = seg == parameter.tasks(task).start.id;
    issf = isosurface(obj, .1);
    k = patch(issf);
    set(k, 'FaceColor', 'r', 'EdgeColor', 'none');
    obj = seg == parameter.tasks(task).end.id;
    issf = isosurface(obj, .1);
    k = patch(issf);
    set(k, 'FaceColor', 'g', 'EdgeColor', 'none');
    for i=1:length(parameter.tasks(task).end.obj.id)
        if parameter.tasks(task).end.obj.id(i) > 0
            obj = seg == parameter.tasks(task).end.obj.id(i);
            issf = isosurface(obj, .1);
            l{i} = patch(issf);
            set(l{i}, 'FaceColor', 'b', 'EdgeColor', 'none');
        end
    end
    view(3);
    daspect([25 25 12]);
    title(['Task: ' num2str(task)]);
    grid on;
    alpha(.5);
    light = camlight('headlight');
    lighting flat;
    legend({'CoM start', 'CoM end', 'CoM split', 'Starting Node & Direction', 'Start Object', 'End Object'});
    pause;
    close all;
end

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
bw3(bw3 ==1) = 255;
for i=1:384
    imwrite(uint8(raw(:,:,i)), [knowledgeDB 'imageSeries/' 'raw' num2str(i, '%.4i') '.png'], 'png');
    imwrite(uint16(seg(:,:,i)), [knowledgeDB 'imageSeries/' 'seg' num2str(i, '%.4i') '.png'], 'png');
    imwrite(uint8(bw3(:,:,i)), [knowledgeDB 'imageSeries/' 'flag' num2str(i, '%.4i') '.png'], 'png');
end

