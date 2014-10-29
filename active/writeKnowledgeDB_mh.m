function [settings, layer, missions] = writeKnowledgeDB_mh(parameter, seg, edges, p, edgesOld, pOld)

% Function to keep all parameter
kdbParameter = setKDBparameter();

% Normalize probability
p = p ./ max(p);

% compute princomps
% objProp = computePrincomps_fun(seg, edges, kdbParameter.settings.voxelSize);
load('/zdata/manuel/sync/activeTraining/20130607pricomps.mat');

% Metadata fuer Tasks generieren
for i=1:size(edges,1) 
    tic;
    % Which edge to look at
    idx1 = edges(i,1);
    idx2 = edges(i,2);  
    % Determine CoM of border between objects (largest area of touch)
    obj1 = seg == idx1;
    obj2 = seg == idx2;
    border = bwdist(obj1) <= 1 & bwdist(obj2) <= 1;
    borderProp = regionprops(border, 'Centroid', 'Area');
    if size(borderProp,1) == 0
        error('No touch detected');
    else
	% Use largest contact area
	[~, idx] = max([borderProp(:).Area]);
	% Local Principal components
	borderCentroid = borderProp(idx).Centroid([2 1 3]).*kdbParameter.settings.voxelSize;
        acceptedPoints1 = find( sum((objProp(idx1).coords - repmat(borderCentroid,[size(objProp(idx1).coords,1) 1])).^2,2) <= 250000);
        [pcbasis1{i},coord{i},weights1{i}] = princomp(objProp(idx1).coords(acceptedPoints1,:));
        sizeObj1{i} = length(acceptedPoints1);
        acceptedPoints2 = find( sum((objProp(idx2).coords - repmat(borderCentroid,[size(objProp(idx2).coords,1) 1])).^2,2) <= 250000);
        [pcbasis2{i},coord2{i},weights2{i}] = princomp(objProp(idx2).coords(acceptedPoints2,:));
        sizeObj2{i} = length(acceptedPoints2);    
        % Isosurface normal
        [gridx, gridy, gridz] = meshgrid([1:size(border,2)], [1:size(border,1)], [1:size(border,3)]);
        gridx = gridx * kdbParameter.settings.voxelSize(1);
        gridy = gridy * kdbParameter.settings.voxelSize(1);
        gridz = gridz * kdbParameter.settings.voxelSize(1);
        iso_border = isosurface(gridx, gridy, gridz, border, 0.5);
        ison_border = isonormals(gridx,gridy,gridz,border,iso_border.vertices);
	% Turn all isonormals to positive
        indneg = find(ison_border(:,1) < 0);
        ison_border(indneg,:) = -ison_border(indneg,:);
	% Turn all isonormals that are on the negative axis
        indneg = find(sum(ison_border == 0 , 2) == 2 & any(ison_border < 0, 2));
        ison_border(indneg,:) = -ison_border(indneg,:);
	% Normalize & normalized mean
        ison_border = ison_border ./ repmat(sqrt(sum(ison_border.^2,2)),[1 3]);
        border_isonormal_local = mean(ison_border,1);
        border_isonormal{i} = border_isonormal_local / norm(border_isonormal_local);        
    end
    t = toc;
    display(['Edge based stuff, Progress: ' num2str(i./size(edges,1)*100, '%.3f') ' %, Time last edge: ' num2str(t) ' seconds']);
end

save('/zdata/manuel/sync/activeTraining20130607temp.mat');

% Load global counter variables
load(kdbParameter.counterLocation);

% Increase counter variables
if restartCounter
	counter.layerID = 1;
	counter.taskID = countMissons; 
else
	counter.layerID = counter.layerID + 1;
	counter.taskID = counter.taskID + countMissions;
end

% Save global counter variables
save(kdbParameter.counterLocation, 'counter');

% Create struct for layer.json
layer.type = 'segmentation';
layer.resolutions = 1;
layer.layerId = counter.layerID;
layer.bbox = paramBG.bboxBig;
layer.class = class(segBig);

% Define folder in which mission & layer json is placed (as well as segmentation in subfolder corresponding to resolution)
rootFolder = [kdbParameter.writeLocation kdbParameter.settings.name '/segmentation/layer' num2str(layer.layerId) '/'];

% Cubes schreiben
writeKnossosRoi([rootFolder '1/'], kdbParameter.settings.name, layer.bbox(:,1)', segBig, class(segBig));

% Write metadata to JSON
savejson('', kdbParameter.settings, [kdbParameter.writeLocation kdbParameter.settings.name '/settings.json']);
savejson('', layer, [rootFolder '/layer.json']);
savejson('', missions, [rootFolder '/missions.json']);

end

