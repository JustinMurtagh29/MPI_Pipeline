function globalCorrSeg(p)

% Fix k=4
files = dir([p.saveFolder, 'correspondences/*04.mat']);
a = vertcat(files(:).name); a = a(:,5:6); for i=1:size(a,1); b(i) = strcmp('04', a(i,:)); end; files = files(b);
corrGlobal1 = [];
corrGlobal2 = [];
for m = 1:length(files)
    load([p.saveFolder, '/correspondences/', files(m).name])     
    load([p.seg.segFolder 'numEl.mat'])
    cubeCoords1 = [str2num(files(m).name(1:2)),str2num(files(m).name(3:4)), str2num(files(m).name(5:6))];
    cubeCoords2 = [str2num(files(m).name(7:8)),str2num(files(m).name(9:10)), str2num(files(m).name(11:12))];
    corrGlobal1 = [corrGlobal1; uint32(result.correspondences(:,1)) + numElTotal(cubeCoords1(1), cubeCoords1(2), cubeCoords1(3))];
    corrGlobal2 = [corrGlobal2; uint32(result.correspondences(:,2)) + numElTotal(cubeCoords2(1), cubeCoords2(2), cubeCoords2(3))];
end

uniqueID = unique([corrGlobal1 corrGlobal2]);
components = findconnectedcomponents1([corrGlobal1 corrGlobal2],uniqueID);
save([p.saveFolder 'seg/m.mat'], 'components');

end
