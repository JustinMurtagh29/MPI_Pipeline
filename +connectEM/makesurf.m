function makesurf(cube, filename)
    issfs = isosurface(cube, .05);
    %issfs = reducepatch(issfs, .01);
    issfs.vertices(:,[1 2]) = issfs.vertices(:,[2 1]);
    issfs.vertices = issfs.vertices .* repmat([11.24 11.24 14],size(issfs.vertices,1),1);
    exportSurfaceToAmira({issfs}, filename);
