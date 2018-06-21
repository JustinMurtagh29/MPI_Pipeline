function dataDir = getDir(dirName)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    thisDir = fileparts(mfilename('fullpath'));
    dataDir = fullfile(thisDir, dirName);
end
