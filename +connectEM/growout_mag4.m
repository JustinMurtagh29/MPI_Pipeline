function [y, mincoord] = growout_mag4 (agglo, segmentMeta)
todos = [0,0,0;0,0,1;0,1,0;0,1,1;1,0,0;1,0,1;1,1,0;1,1,1];
base_c = segmentMeta.point(agglo,:) / 128 / 4;
cube_coord = [];
for idx = 1 : 8
    cube_coord = [cube_coord; floor(base_c) + bsxfun(@times, 1 - 2 * ((base_c - floor(base_c)) < 0.25), todos(idx, :))];
end
y = false((max(cube_coord) - min(cube_coord) + 1) * 128);
cube_coord_unique = unique(cube_coord, 'rows');
agglo_double =  repmat(agglo, 8, 1);
for idx = 1 : size(cube_coord_unique, 1)
    idx
    thiscube = readKnossosCube('/gaba/wKcubes/Connectomics department/2012-09-28_ex145_07x2_ROI2017/segmentation/4/', ...
    '2012-09-28_ex145_07x2_ROI2016_corrected_mag4', cube_coord_unique(idx, :), 'uint32');

    hits = find(ismember(cube_coord, cube_coord_unique(idx, :), 'rows'));
    a = @(n)cube_coord_unique(idx, n) - min(cube_coord_unique(:, n));
    b = @(n)a(n) * 128 + 1 : (a(n) + 1) * 128;
    y(b(1), b(2), b(3)) = ismember(thiscube, agglo_double(hits));
    

%        if sum(sum(sum(thiscube == agglo(hits(idx_h))))) == 0
%            warning(num2str(agglo(hits(idx_h))));
%        else
%            connectEM.makesurf(ytemp, [num2str(agglo(hits(idx_h))), '.issf']);
%        end
%    end
end
mincoord = min(cube_coord_unique)*128*4;