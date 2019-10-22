function numericToHdf5(file, dset, data)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    assert(isnumeric(data));
    sz = size(data);
    
    % NOTE(amotta): Remove trailing singleton dimension
    if numel(sz) == 2 && sz(2) == 1; sz = sz(1); end
    
    % NOTE(amotta): Compression is counter-productive in this particular
    % case. Storage usage reported by h5ls is around 30 % to 50 % with
    % deflate (9).
    h5create(file, dset, sz, 'Datatype', class(data));
    h5write(file, dset, data);
end
