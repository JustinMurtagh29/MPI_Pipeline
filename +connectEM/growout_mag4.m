function y = growout_mag4 (agglo, segmentMeta)
cube_coord = floor(segmentMeta.point(:, agglo)' / 128 / 4);
y = zeros((max(cube_coord) - min(cube_coord) + 3) * 128);
cube_coord_unique = unique(cube_coord, 'rows');

for idx = 1 : size(cube_coord_unique, 1)
    thiscube = readKnossosRoi('/gaba/wKcubes/Connectomics department/2012-09-28_ex145_07x2_ROI2017/segmentation/4/', ...
    '2012-09-28_ex145_07x2_ROI2016_corrected_mag4', [(cube_coord_unique(idx, :)' - 1) * 128 + 1, (cube_coord_unique(idx, :)' + 2) * 128], 'uint32');
    hits = find(ismember(cube_coord, cube_coord_unique(idx, :), 'rows'));
    for idx_h = 1 : length(hits)
        a = @(n)cube_coord(hits(idx_h), n) - min(cube_coord(:, n)) + 1;
        b = @(n)(a(n)-1) * 128 + 1 : (a(n) + 2) * 128;
        c = {b(1), b(2), b(3)};
        ytemp = y;
        ytemp = ytemp * 0;
        ytemp(c{:}) = ytemp(c{:}) | thiscube == agglo(hits(idx_h));

        y(c{:}) = y(c{:}) | thiscube == agglo(hits(idx_h));


        if sum(sum(sum(thiscube == agglo(hits(idx_h))))) == 0
            warning(num2str(agglo(hits(idx_h))));
%        else
%            connectEM.makesurf(ytemp, [num2str(agglo(hits(idx_h))), '.issf']);
        end
    end
end
