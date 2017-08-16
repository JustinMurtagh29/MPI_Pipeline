baseFolder = '/tmpscratch/kboerg/';
prefix = 'visX24_';

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

save('/tmpscratch/mberning/superagglos.mat', 'superagglos');

% Also save as table
fH = @(x)array2table(x, 'VariableNames', {'x' 'y' 'z' 'segmentID'});
for i=1:length(superagglos)
    superagglos(i).nodes = fH(superagglos(i).nodes);
end
save('/tmpscratch/mberning/superagglosTable.mat', 'superagglos');

clear temp;
% And as graph
for i=1:length(superagglos)
    temp{i} = graph(table(superagglos(i).edges, 'VariableNames', {'EndNodes'}), superagglos(i).nodes);
end
superagglos = temp;
clear temp;
save('/tmpscratch/mberning/superagglosGraph.mat', 'superagglos');

