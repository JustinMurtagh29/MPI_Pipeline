function arrayToHdf5(file, dsetOrGroup, array)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    if isnumeric(array)
        numericToHdf5(file, dsetOrGroup, array);
    elseif iscell(array)
        cellToHdf5(file, dsetOrGroup, array);
    elseif iscategorical(array)
        categoricalToHdf5(file, dsetOrGroup, array)
    elseif isstruct(array)
        structToHdf5(file, dsetOrGroup, array);
    else
        error('Unhandled class "%s"', class(array));
    end
end
