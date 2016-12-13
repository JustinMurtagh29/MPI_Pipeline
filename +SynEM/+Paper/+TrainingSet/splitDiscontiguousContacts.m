function [ contactSurfaceList, neighborIDsNew ] = splitDiscontiguousContacts( contacts,neighborIDs, siz,bbox,useInGUI )
%SPLITDISCONTIGUOUSCONTACTS Split discontigous contacts determined by
%regionprops.

fprintf('%s INFO splitting discontiguous contacts\n',datestr(now,'yyyy-mm-dd HH:MM:SS,FFF'));

%preallocate arrays for speed
contactSurfaceList = cell(4*length(contacts),1);
neighborIDsNew = zeros(4*length(contacts),2);

contactCount = 1;
for i=1:length(contacts)
    %only consider bounding box of current contact
    [x,y,z] = ind2sub(siz,contacts{i});
    u = bsxfun(@minus,x,min(x(:)) - 1);
    v = bsxfun(@minus,y,min(y(:)) - 1);
    w = bsxfun(@minus,z,min(z(:)) - 1);
    U = [u,v,w];
    contactSurface = false(max(u),max(v),max(w));
    indices = sub2ind(size(contactSurface),U(:,1),U(:,2),U(:,3));
    contactSurface(indices) = true;
    L = bwlabeln(contactSurface);
    STATS = regionprops(L,'PixelIdxList');
    for k=1:length(STATS)
        if length(STATS(k).PixelIdxList) > 99
            indices = STATS(k).PixelIdxList;
            [r,s,t] = ind2sub(size(contactSurface),indices);
            r = bsxfun(@plus,r,min(x(:)) - 1);
            s = bsxfun(@plus,s,min(y(:)) - 1);
            t = bsxfun(@plus,t,min(z(:)) - 1);
            com = mean([s,r,t]);
            wallIndices = [s,r,t];
            idx = knnsearch(wallIndices,com);
            centroid = wallIndices(idx,:);
            centroid = [centroid(2),centroid(1),centroid(3)];
            %record contact if in inner bounding box
            if isInInnerBox(centroid,bbox,useInGUI);
                contactSurfaceList{contactCount} = sub2ind(siz,r,s,t);
                neighborIDsNew(contactCount,:) = neighborIDs(i,:);
                contactCount = contactCount + 1;
            end
        end
    end
end

contactSurfaceList(contactCount:end) = [];
neighborIDsNew(contactCount:end,:) = [];

end

function isValidContactSurface = isInInnerBox( centerOfMass,bbox,useInGUI )
%CHECKINNERBOX Check of centroid lies in a box which is smaller by
%[100,100,40] at each boundary of the bbox
if useInGUI
    innerBox = [100,bbox(1,2)-bbox(1,1)-99;100,bbox(2,2)-bbox(2,1)-99;40,bbox(3,2)-bbox(3,1)-39];
else
    innerBox = [40,bbox(1,2)-bbox(1,1)-39;40,bbox(2,2)-bbox(2,1)-39;14,bbox(3,2)-bbox(3,1)-13]; %[24,24,8] voxels for the filter boundaries and [14,14,5] voxels for the centroid to be not too close to the boundary
end
isValidContactSurface = (innerBox(1,1) < centerOfMass(1)) & (centerOfMass(1) < innerBox(1,2)) & ...
                (innerBox(2,1) < centerOfMass(2)) & (centerOfMass(2) < innerBox(2,2)) & ...
                (innerBox(3,1) < centerOfMass(3)) & (centerOfMass(3) < innerBox(3,2));

end

