oldDir = '/zdata/manuel/results/parameterSearch/';
newDir = '/zdata/manuel/parameterSearchTemp/';

files = dir(oldDir);
files(1:2) = []; % remove . and ..
for i=1:length(files)
    mkdir([newDir files(i).name '/']);
    filesSubfolder = dir([oldDir files(i).name '/iter*']);
    for j=1:length(filesSubfolder)
        mkdir([newDir files(i).name '/' filesSubfolder(j).name '/']);
        filesSubSubfolder = dir([oldDir files(i).name '/' filesSubfolder(j).name '/gpu*']);
        for k=1:length(filesSubSubfolder)
            matFiles = dir([oldDir files(i).name '/' filesSubfolder(j).name '/' filesSubSubfolder(k).name '/*.mat']);
            if ~isempty(matFiles)
                mkdir([newDir files(i).name '/' filesSubfolder(j).name '/' filesSubSubfolder(k).name '/']);
                copyfile([oldDir files(i).name '/' filesSubfolder(j).name '/' filesSubSubfolder(k).name '/' matFiles(end).name], [newDir files(i).name '/' filesSubfolder(j).name '/' filesSubSubfolder(k).name '/' matFiles(end).name]);
            end
        end
    end
end

