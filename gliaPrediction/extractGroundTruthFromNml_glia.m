function [labeledIdx, labels] = extractGroundTruthFromNml(seg, segments, skelParam)

skel = skeleton(skelParam.file,0);
bbox = skelParam.bbox;

%% get tree ids

gliaTrees = [];
for i=1:length(skel.nodesAsStruct)
    for j=1:length(skel.nodesAsStruct{i})
        if strcmp(skel.nodesAsStruct{i}(j).comment,'glia') || strcmp(skel.nodesAsStruct{i}(j).comment,'Glia') 
            gliaTrees(end+1) = skel.thingIDs(i);
            continue;
        end
    end
end

idxBbox = find(cellfun(@(x) strcmp('bbox',x),skel.names));

%% get glia and nonGlia ids
glia.ids = [];
nonGlia.ids = [];

for i=1:length(skel.nodes)
    if i == idxBbox
        continue;
    end
    gliaFlag = any(gliaTrees == skel.thingIDs(i));
    for j=1:size(skel.nodes{i})
        %restrict on bbox, switch to local coords
        if(any(bbox(:,1)' - skel.nodes{i}(j,1:3) > 0) || any(bbox(:,2)' - skel.nodes{i}(j,1:3) < 0))
            continue;
        end
        localCoords = transformCoords(skel.nodes{i}(j,1:3),bbox,0);
        if gliaFlag
            glia.ids(end+1) = seg(localCoords(1),localCoords(2),localCoords(3));
        else
            nonGlia.ids(end+1) = seg(localCoords(1),localCoords(2),localCoords(3));
        end           
    end
end

% resolve id ambiguities if ratio is small enough
commonIds = intersect(glia.ids,nonGlia.ids);
for i=1:length(commonIds)

    N(i) = sum(glia.ids == commonIds(i));
    G(i) = sum(nonGlia.ids == commonIds(i));
    ratio(i) = min(N(i),G(i))/max(N(i),G(i));
end
    nodesPerSeg = [N;G;ratio;commonIds];  
    corrected = nodesPerSeg(:,nodesPerSeg(3,:) < 0.15);

for i=1:length(corrected)
    if corrected(1,i) > corrected(2,i)
        nonGlia.ids(nonGlia.ids == corrected(4,i)) = [];
    else
        glia.ids(glia.ids == corrected(4,i)) = [];
    end
end  

glia.ids = unique(glia.ids);
nonGlia.ids = unique(nonGlia.ids);

commonIds = intersect(glia.ids,nonGlia.ids);
glia.ids = setdiff(glia.ids,commonIds);
nonGlia.ids = setdiff(nonGlia.ids,commonIds);
% exclude 0
nonGlia.ids(nonGlia.ids == 0) = [];


%% Find labels

segmentIds = cell2mat({segments.id});

for i = 1:length(segmentIds)
    if any(glia.ids == segmentIds(i))
        labels(i) = 1;
    elseif any(nonGlia.ids == segmentIds(i))
        labels(i) = -1;
    else
        labels(i) = 0;
    end
end

labeledIdx = labels ~= 0;
labels(labels == 0) = [];

% Transpose in order to have same dimension as supervoxel weights
labeledIdx = labeledIdx';
labels = labels';

end

