function princomps = computePrincomps_fun(segNew, edges, resolution)

objlist = unique(edges);
princomps = struct('objId', {}, 'objSize', {}, 'objProp', {}, 'pcBase1', {}, 'pcBase2', {}, 'pcBase3', {}, 'pcSpan1', {}, 'pcSpan2', {}, 'pcSpan3', {});

for i=objlist'
    tic;
    princomps(i).objectId = i;
    objProp = regionprops(segNew == i, 'PixelList', 'PixelIdxList', 'Area', 'Centroid');
    princomps(i).objSize = objProp.Area;
    princomps(i).objPos = objProp.Centroid.*resolution;
    % coordinates in nm
    princomps(i).coords = objProp.PixelList.*repmat(resolution, [size(objProp.PixelList,1) 1]);     
    % compute princomp
    [pcbasis,coord_in_new_basis,weights] = princomp(princomps(i).coords);
    % add to list
    princomps(i).pcBase1 = pcbasis(:,1);
    princomps(i).pcBase2 = pcbasis(:,2);
    princomps(i).pcBase3 = pcbasis(:,3);
    princomps(i).pcSpan1 = max(coord_in_new_basis(:,1)) - min(coord_in_new_basis(:,1));
    princomps(i).pcSpan2 = max(coord_in_new_basis(:,2)) - min(coord_in_new_basis(:,2));
    princomps(i).pcSpan3 = max(coord_in_new_basis(:,3)) - min(coord_in_new_basis(:,3));
    t = toc;
    display(['Object based stuff, Progress: ' num2str(i ./ max(objlist) * 100, '%.3f') ' %, Time last edge: ' num2str(t) ' seconds']);
end

end

