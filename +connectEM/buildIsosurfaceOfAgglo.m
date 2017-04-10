function buildIsosurfaceOfAgglo(p, segIds, fileName)

    if ~isempty(segIds)
        issfs = Visualization.buildIsoSurface(p, segIds, 'reduce', 0.05);
        issfs.vertices = bsxfun(@times, issfs.vertices, p.raw.voxelSize);
        Visualization.writePLY(issfs, [1, 0, 0], fileName);
    end

end

