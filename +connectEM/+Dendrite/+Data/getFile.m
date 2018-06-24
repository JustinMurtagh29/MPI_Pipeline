function path = getFile(name)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    thisDir = fileparts(mfilename('fullpath'));
    path = fullfile(thisDir, name);
end
