function collectSuperagglos(baseFolder,prefix,outputFolder)
if ~exist('baseFolder','var')
    baseFolder = '/tmpscratch/kboerg/';
end
if ~exist('prefix','var')
    prefix = 'visX24_';
end
if ~exist('outputFolder','var')
    outputFolder = '/tmpscratch/mberning/';
end
superagglos = [];
folders1 = dir(fullfile(baseFolder, strcat(prefix, '*')));
for i=1:length(folders1)
    folders2 = dir(fullfile(baseFolder, folders1(i).name, strcat(prefix, '*')));
    for j=1:length(folders2)
        idx_agglo = str2double(folders2(j).name(8:end));
        temp = load(fullfile(baseFolder, folders1(i).name, folders2(j).name, 'superagglo.mat'));
        superagglos(idx_agglo).nodes = cat(2, temp.nodes2, temp.segIds2);
        superagglos(idx_agglo).edges = unique(sort(temp.edges2,2), 'rows');
    end
end

save(fullfile(outputFolder,'superagglos.mat'), 'superagglos');

% Also save as table
fH = @(x)array2table(x, 'VariableNames', {'x' 'y' 'z' 'segmentID'});
for i=1:length(superagglos)
    superagglos(i).nodes = fH(superagglos(i).nodes);
end
save(fullfile(outputFolder,'superagglosTable.mat'), 'superagglos');

clear temp;
% And as graph
for i=1:length(superagglos)
    temp{i} = graph(table(superagglos(i).edges, 'VariableNames', {'EndNodes'}), superagglos(i).nodes);
end
superagglos = temp;
clear temp;
save(fullfile(outputFolder,'superagglosGraph.mat'), 'superagglos');

