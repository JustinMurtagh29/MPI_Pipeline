load('/gaba/u/mberning/results/pipeline/20151111T183414/allParameter.mat');
load('/gaba/u/mberning/2012-09-28_ex145_07x2_skeletons/skelGT.mat');
graph = load([p.saveFolder 'graph.mat']);

prob = [1:-0.01:.91];
excludeIds = cat(1,gt.segIds{:});
for p_i=1:length(prob)
    display(['Agglomeration: ' num2str(prob(p_i), '%3.2f')]);
    tic;
    for i=1:length(gt.segIds)
        if p_i == 1
            ids{i,p_i} = gt.segIds{i};
            nrAgglo(i,p_i) = 0;
        else
            [ids{i,p_i}, nrAgglo(i,p_i)] = agglomerateSG_simple(graph, ids{i, p_i-1}, prob(p_i), excludeIds);
        end
        Util.progressBar(i,length(gt.segIds));
    end
    display(['Segments agglomerated: ' num2str(sum(nrAgglo(:,p_i)))]);
    toc;
    % Check for & remove overlaps if any are found
    display('Checking for overlaps');
    tic;
    overlap = cell(length(gt.segIds));
    for k=1:length(gt.segIds)
        for j=k+1:length(gt.segIds)
            overlap{k,j} = intersect(ids{k,p_i},ids{j,p_i});
        end
    end
    newExcludeIds = unique(cat(1,overlap{:}));
    excludeIds = cat(1,excludeIds,newExcludeIds);
    keepIdx = cellfun(@(x)~ismember(x, newExcludeIds), ids(:,p_i),'UniformOutput', false);
    ids(:,p_i) = cellfun(@(x,y)x(y), ids(:,p_i), keepIdx, 'UniformOutput', false);
    nrExcluded = sum(cellfun(@(x)sum(~x),keepIdx));
    display(['Segments removed due to overlap: ' num2str(nrExcluded) ', unique: ' num2str(length(newExcludeIds))]);
    toc;
end
clear i p_i;
save([p.saveFolder 'agglomeration/' datestr(clock,30) '.mat']);

%% Write isosurfaces for visualization in Amira/Blender

% First resort to have better numbering for Heiko
[gtSorted.label, idxResort] = sort(gt.label);
gtSorted.segIds = gt.segIds(idxResort);
idsSorted = ids(idxResort,:);

p_i = 4; % Chose 97% agglomeration for now
for i=1:length(gtSorted.label)
    inputCell{i} = {p, idsSorted{i,p_i}, ...
        ['/gaba/u/mberning/2012-09-28_ex145_07x2_skeletons/isoFull/new/' num2str(prob(p_i), '%3.2f') ...
        '/' gtSorted.label{i} '_' num2str(i, '%.5i')  '.issf'], ...
        ['/gaba/u/mberning/2012-09-28_ex145_07x2_skeletons/isoFull/new/' num2str(prob(p_i), '%3.2f') ...
        '/lowRes_' gtSorted.label{i} '_' num2str(i, '%.5i')  '.issf']};
end
clear i p_i;
functionH = @galleryCortexCCfromSG;

CLUSTER_CPU.IndependentSubmitFcn{2} = '-p -400 -pe openmp 1 -l h_vmem=32G,h_rt=100:00:00';
startCPU(functionH, inputCell, ['agglomeration isosurfaces for Amira']);
clear functionH inputCell;


% Quick and dirty subsampling part two
for r=1
    idx = 1;
    for p_i=5
        for i=1:length(gt.segIds)
            inputFile = ['/gaba/u/mberning/2012-09-28_ex145_07x2_skeletons/isoFull/' num2str(prob(p_i), '%3.2f') ...
                '/' resolutions{r} '_' gt.label{i} '_' num2str(i, '%.5i')  '.mat']; 
            outputFile = ['/gaba/u/mberning/2012-09-28_ex145_07x2_skeletons/isoFull/' num2str(prob(p_i), '%3.2f') ...
                '_smaller/' resolutions{r} '_' gt.label{i} '_' num2str(i, '%.5i')  '.issf'];
            inputCell{r}{idx} = {inputFile outputFile};
            idx = idx + 1;
        end
    end
end
functionH = @subsampleAgain;

