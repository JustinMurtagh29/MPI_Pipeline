load('/gaba/u/mberning/2012-09-28_ex145_07x2_skeletons/skelGT.mat');

prob = [1:-0.01:.96];
excludeIds = cat(1,gt.segIds{:});
% Exclude vessels for now, could take long
gt.segIds(end) = [];
gt.label(end) = [];
for p_i=6:length(prob)
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
    ids(:,p_i) = cellfun(@(x)x(~ismember(x, newExcludeIds)), ids(:,p_i),'UniformOutput', false);
    display(['Segments removed due to overlap: ' num2str(length(newExcludeIds))]);
    toc;
end

resolutions = {'isoSurfLowRes' 'isoSurfMidRes' 'isoSurf'};
for r=1:length(resolutions)
    idx = 1;
    for p_i=1:length(prob)
        for i=1:length(gt.segIds)
            inputCell{r}{idx} = {p, ids{i,p_i}, ...
                ['/gaba/u/mberning/2012-09-28_ex145_07x2_skeletons/isoFull/' num2str(prob(p_i), '%3.2f') ...
                '/' resolutions{r} '_' gt.label{i} '_' num2str(i, '%.5i')  '.issf'], resolutions{r}};
            idx = idx + 1;
        end
    end
end
functionH = @galleryCortexCCfromSG_new;

CLUSTER_CPU.IndependentSubmitFcn{2} = '-p -400 -pe openmp 1 -l h_vmem=16G,h_rt=100:00:00';
for r=1:length(resolutions)
    startCPU(functionH, inputCell{r}, ['agglomeration full isosurfaces' resolutions{r}]);
end
clear p_i i idx functionH inputCell;

