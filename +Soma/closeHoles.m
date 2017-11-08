function [ newSegIds, name ] = closeHoles( p, meta, rp, graph, voxelSize, somaID, aggloSegIds )
%CLOSEHOLES Close holes in an soma agglomerate by searching for surface normals 
%   that have an opposite direction to the vector from the center of the
%   nuclei to the point of interest.
   
    nrSlices = 8;

    %get agglo bbox
    boxes = Seg.Global.getSegToBoxMap(p); 
    bboxes = boxes(:,:,aggloSegIds);
    xMinG = min(bboxes(1,1,:));
    yMinG = min(bboxes(2,1,:));
    zMinG = min(bboxes(3,1,:));
    xMaxG = max(bboxes(1,2,:));
    yMaxG = max(bboxes(2,2,:));
    zMaxG = max(bboxes(3,2,:));
    sliceSize = ceil((zMaxG-zMinG)/nrSlices);
    
    %get center point coordinates
    mag1bbox = [128, 128, 128, 5573, 8508, 3413]; mag1bbox = Util.convertWebknossosToMatlabBbox(mag1bbox);
    mag1bbox = mag1bbox + [-25 25; -25 25; -10 10];
    centerPoint = bsxfun(@plus,4.* rp(somaID).Centroid, mag1bbox(:,1)');
    
    newSegIds = [];
    
    %every slice is cat into 4 quarters to have smaller bounding boxes
    for j=0:nrSlices-1
        disp(j)
        for q = 1:4
            if q == 1
                xMin = xMinG;
                yMin = yMinG;
                zMin = zMinG;
                xMax = xMinG + ceil((xMaxG - xMinG)/2);
                yMax = yMinG + ceil((yMaxG - yMinG)/2);
            elseif q == 2
                xMin = xMinG;
                yMin = yMinG + ceil((yMaxG - yMinG)/2);
                zMin = zMinG;
                xMax = xMinG + ceil((xMaxG - xMinG)/2);
                yMax = yMax;
            elseif q == 3
                xMin = xMinG + ceil((xMaxG - xMinG)/2);
                yMin = yMinG;
                zMin = zMinG;
                xMax = xMaxG;
                yMax = yMinG + ceil((yMaxG - yMinG)/2);
            elseif q == 4
                xMin = xMinG + ceil((xMaxG - xMinG)/2);
                yMin = yMinG + ceil((yMaxG - yMinG)/2);
                zMin = zMinG;
                xMax = xMaxG;
                yMax = yMaxG;
            end

            %get segmentation with new bounding boxes
            somaSeg = readKnossosRoi('/gaba/wKcubes.new/Connectomics department/2012-09-28_ex145_07x2_ROI2017/segmentation/1','2012-09-28_ex145_07x2_ROI2016_corrected_mag1', ...
                [xMin-1 xMax+1; yMin-1 yMax+1; ((zMin-1)+sliceSize*j) ((zMin-1)+sliceSize*(j+1))], 'uint32', '', 'raw');

            %get mask for inverse soma. Somehow normals seem to be correctly when not using the soma but the inverse. Should be fixed.
            soma = ~ismember(somaSeg, aggloSegIds);
            %this line allocates a lot of ram and is the reason for the slices
            soma = Seg.Local.closeGaps(soma);
            disp('critical part survived');

            soma = smooth3(soma,'box',9);
            [m,n,o] = size(soma);
            [X,Y,Z] = meshgrid(1:n,1:m,1:o);
            somaIso = isosurface(X,Y,Z,soma,.5,'noshare');
            disp('end iso')
            
            isoVertices = somaIso.vertices;
            %check if isoVertices is empty
            if size(isoVertices, 1)==0
                continue;
            end
            
            %remove borders from cutting
            %lowest values are 1.5
            isoVertices(any(isoVertices==1.5,2),:) = [];
            if size(isoVertices, 1)==0
                continue;
            end
            %highest are calculated
            xMax2 = xMax+1-xMin-1+0.5;
            yMax2 = yMax+1-yMin-1+0.5;
            zMax2 = (zMin-1)+sliceSize*(j+1)-(zMin-1)+sliceSize*j+0.5;
                        size(isoVertices)

            isoVertices(isoVertices(:,1)==yMax2,:) = [];
            if size(isoVertices, 1)==0
                continue;
            end
            isoVertices(isoVertices(:,2)==xMax2,:) = [];
             if size(isoVertices, 1)==0
                continue;
            end
            isoVertices(isoVertices(:,3)==zMax2,:) = [];
             if size(isoVertices, 1)==0
                continue;
            end

            somaIso.vertices = isoVertices;

            %randomly take 50000 vertices
            if size(somaIso.vertices,1) < 50000
                nrSamples = size(somaIso.vertices,1);
            else
                nrSamples = 50000;
            end
            downSampleIdx = randperm(size(somaIso.vertices,1),nrSamples);
            vertices = somaIso.vertices(downSampleIdx,:);

            %calculate normals
            isoNormals = isonormals(soma,vertices);
            isoNormals = bsxfun(@times,isoNormals,voxelSize);
            isoNormals = normr(isoNormals);

            %calculate angle between normal of point of interest and centroid of soma
            angles = zeros(size(downSampleIdx,2),1);
            for i=1:size(downSampleIdx,2)
                disp(i)
                refPoint = vertices(i,:);
                refPoint([1,2]) = refPoint([2,1]);
                refPoint = refPoint + [xMin, yMin, zMin];
                refVec = refPoint-centerPoint;
                refVec = refVec .* voxelSize;
                refVec([1,2]) = refVec([2,1]);
                refVec = refVec .* -1;
                angles(i) = atan2(norm(cross(isoNormals(i,:)',refVec')),dot(isoNormals(i,:)',refVec'));    
            end

            %hold off
            %scatter3(somaIso.vertices(:,1),somaIso.vertices(:,2),somaIso.vertices(:,3))
            %hold on

            %keep all Vertices and normals that build an angle > x
            activeVertices = vertices(angles>2.9,:);
            activeNormals = isoNormals(angles>2.9,:);

            %scatter3( vertices(angles>2.9,1),  vertices(angles>2.9,2),  vertices(angles>2.9,3));
            %scatter3(soma(downSampleIdx(find(angles>2.4)),1),soma(downSampleIdx(find(angles>2.4)),2),soma(downSampleIdx(find(angles>2.4)),3));

            %change x and y because they are changed by isosurface
            activeVertices(:,[1,2]) = activeVertices(:,[2,1]);
            activeNormals(:,[1,2]) = activeNormals(:,[2,1]);

	    %randomly take 3000 active vertices
            if size(activeVertices,1) < 3000
                nrSamples2 = size(activeVertices,1);
            else
                nrSamples2 = 3000;
            end
            downSampleIdx = randperm(size(activeVertices,1),nrSamples2);
            %downSampleIdx = randperm(size(somaIso.vertices,1),5000);
            activeVertices = activeVertices(downSampleIdx,:);
            activeNormals = activeNormals(downSampleIdx,:);

            %get new seg ids
            %search length determines how far the picking up of segments reaches. When along the search path a somatic segment is found all segments are picked up. If not all are removed. This is done because dentrites potentially also build opposing angles but have a longer distance to the soma.
            searchLength = 20;
            addedSegIds = [];
            for i=1:size(activeVertices,1)
                disp(i);
                actualVertice = activeVertices(i,:);
                actualNormal = activeNormals(i,:)*-1;
                segIds = zeros(searchLength,1);
                for s=1:searchLength
                    try
                        coord = round(actualVertice + s*actualNormal);
                        segIds(s) = somaSeg(coord(1), coord(2), coord(3));
                    end
                end
                segIds = unique(segIds);
                segIds(segIds==0)=[];
                if size(intersect(segIds, aggloSegIds),1) > 0 & size(setdiff(segIds, aggloSegIds)) > 0
                    addedSegIds = cat(1, addedSegIds, setdiff(segIds, aggloSegIds));
                end
            end
            addedSegIds = unique(addedSegIds);
            newSegIds = cat(1, newSegIds, addedSegIds);
            newSegIds = unique(newSegIds);
        end
    end
    %fill all holes that are completely surrounded by agglomerate
    newSegIds = cat(1, aggloSegIds, newSegIds);
    newSegIds = unique(newSegIds);
    surridx = any(ismember(graph.edges, newSegIds),2);
    theseEdges = graph.edges(~surridx,:);
    CCcut = Graph.findConnectedComponents(theseEdges, false, true);
    [maxsize] = cellfun(@(C) size(C,1), CCcut);
    CCcut(maxsize==max(maxsize)) = [];
    if size(CCcut) ~= 0
        for i=1:size(CCcut)
            newSegIds=cat(1, newSegIds, CCcut{i});
        end
    end
    newSegIds=unique(newSegIds);
    name = strcat('soma',int2str(somaID));
end

