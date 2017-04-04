function buildIsosurafceOfAgglo(p, segIds, fileName)

    issfs = Visualization.buildIsoSurface(p, segIds, 'reduce', 0.05);
    issfs{1}.vertices = bsxfun(@times, issfs{1}.vertices, p.raw.voxelSize);
    Visualization.writePLY(issfs{1}, fileName{1});

end

