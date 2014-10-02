function globalCorrSeg(p)
%find corresponding segments for all cubes w/ connComp
% Fix k=4
files = dir([p.correspondence.saveFolder, '*.mat']);
corrNewCoords1 = [];
corrNewCoords2 = [];
for m = 1:length(files)
    
    load([p.correspondence.saveFolder files(m).name])     
    load([p.seg.root 'numEl.mat'])
    cubeCoords1 = [str2num(files(m).name(1:2)),str2num(files(m).name(3:4)), str2num(files(m).name(5:6))];
    cubeCoords2 = [str2num(files(m).name(7:8)),str2num(files(m).name(9:10)), str2num(files(m).name(11:12))];
    load([p.local(cubeCoords1(1), cubeCoords1(2), cubeCoords1(3)).saveFolder, 'globalSegId.mat']);
    newSeg1 = checkMat;
    load([p.local(cubeCoords2(1), cubeCoords2(2), cubeCoords2(3)).saveFolder, 'globalSegId.mat']);
    newSeg2 = checkMat;

    %calculate globalIDs for correspondences
    result.correspondences(:,1) = uint32(result.correspondences(:,1)) + numElTotal(cubeCoords1(1), cubeCoords1(2), cubeCoords1(3));
    result.correspondences(:,2) = uint32(result.correspondences(:,2)) + numElTotal(cubeCoords2(1), cubeCoords2(2), cubeCoords2(3));
    
    for i = 1:size(result.correspondences,1)
        idx1 = find(newSeg1(:,1) == result.correspondences(i,1));
        [x,y] = ind2sub([size(newSeg1,1), size(newSeg1,2)], idx1);
   
        corrNewCoords1 = [corrNewCoords1; newSeg1(x,y)];
        idx2 = find(newSeg2(:,1) == result.correspondences(i,2));
        [x,y] = ind2sub([size(newSeg2,1), size(newSeg2,2)], idx2);
        corrNewCoords2 = [corrNewCoords2; newSeg2(x,y)];
    end
end

corrGlobal = [corrNewCoords1, corrNewCoords2];
uniqueID = unique([corrGlobal]);
components = findconnectedcomponents1([corrNewCoords1 corrNewCoords2],uniqueID);
save([p.seg.root 'globalMapping.mat'], 'components');

%create mapping vector for Oxalis
maxElements = numElTotalAll(end,end,end);
mapVec = unit32(1:maxElements);
mapVec(corrNewCoords1)  = corrNewCoords2;
save([p.seg.root 'correspondenceMapping.mat'], 'mapVec');

end
