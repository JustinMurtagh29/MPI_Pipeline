function globalCorrSeg(p)

files = dir([p.saveFolder, 'correspondences/*.mat']);

    for m = 1:length(files)

        load([p.saveFolder, '/correspondences/', files(m).name])     
        load([p.seg.segFolder 'numEl.mat'])
        cubeCoords1 = [str2num(files(m).name(2)),str2num(files(m).name(4)), str2num(files(m).name(6))]
        cubeCoords2 = [str2num(files(m).name(8)),str2num(files(m).name(10)), str2num(files(m).name(12))];

        corrGlobal(:,1) = result.correspondences(:,1) + numEl(cubeCoords1(1), cubeCoords1(2), cubeCoords1(3));
        corrGlobal(:,2) = result.correspondences(:,2) + numEl(cubeCoords2(1), cubeCoords2(2), cubeCoords2(3));

    end
end
