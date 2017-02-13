function globalCorrSeg(p, fname)

    % Not very pretty, parsing involved cubes from filenames xD
    cubeCoords1 = [str2num(fname(1:2)),str2num(fname(3:4)), str2num(fname(5:6))];
    cubeCoords2 = [str2num(fname(7:8)),str2num(fname(9:10)), str2num(fname(11:12))];
    % Load correspondences
    load([p.correspondence.saveFolder fname]);
    % Globalize
    load([p.local(cubeCoords1(1),cubeCoords1(2),cubeCoords1(3)).saveFolder 'localToGlobalSegId.mat']);
    corrLocal1 = uint32(uniqueCorrespondences(:,1));
    corrGlobal1 = zeros(size(corrLocal1),'uint32');
    for i=1:length(corrLocal1)
        corrGlobal1(i) = globalIds(localIds == corrLocal1(i));
    end 
    load([p.local(cubeCoords2(1),cubeCoords2(2),cubeCoords2(3)).saveFolder 'localToGlobalSegId.mat']);
    corrLocal2 = uint32(uniqueCorrespondences(:,2));
    corrGlobal2 = zeros(size(corrLocal2),'uint32');
    for i=1:length(corrLocal2)
        corrGlobal2(i) = globalIds(localIds == corrLocal2(i));
    end
    % Save globalized version alongside old one
    Util.save([p.correspondence.saveFolder fname(1:end-4) 'global.mat'], corrGlobal1, corrGlobal2);

end
