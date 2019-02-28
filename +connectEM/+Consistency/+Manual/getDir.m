function dataDir = getDir(dirName)
    thisDir = fileparts(mfilename('fullpath'));
    dataDir = fullfile(thisDir, dirName);
end
