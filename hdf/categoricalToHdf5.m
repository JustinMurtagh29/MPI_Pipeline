function categoricalToHdf5(outFile, dset, cats)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    assert(iscategorical(cats));
    catNames = categories(cats);
    
    assert(numel(catNames) <= intmax('uint8'));
    numericToHdf5(outFile, dset, uint8(cats));
    
    for curIdx = 1:numel(catNames)
        curName = lower(catNames{curIdx});
        h5writeatt(outFile, dset, curName, uint8(curIdx));
    end
end
