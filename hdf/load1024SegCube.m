function data = load1024SegCube(p)
    % point close to or at center of dataset
    % taken from webKNOSSOS (segNew)
    box = [5121, 3584, 1728];
    box = [box(:) - 512, box(:) + 511];

    data = loadSegDataGlobal(p.seg.root, p.seg.prefix, box);
end