function addQueries(agglos,outputFolder,skeletonFolders,segmentsLeftover,isaxon)
% this function applies the queries on the agglos/superagglos (new
% representation!)
%
% INPUT
% agglos            agglomerations in the new representation (nodes/edges)
% outputFolder      what it says
% skeletonFolders   cell array of strings with the absolute path to folders
%                   containing the queries in nml format
% segmentLeftover   if agglo does only contain a subpart of all agglos
%                   (e.g. above 5 micron length), the other can be added here
% isaxon            boolean if the agglos are from the axon class
%                   (additional steps done there)



%%%%%
% outputFolder = '/tmpscratch/mberning/axonQueryResults/';
% outputFolder ='/gaba/scratch/kboerg/queryResultsDendrite20170717/';
% load('/gaba/scratch/mberning/axonQueryGeneration/beforeQueryGeneration.mat', 'axonsNew');
%%%%%
%%% old axon round
% scratchFolder = outputFolder;
% skeletonFolders = {'MBKMB_L4_axons_queries_2017_a' 'MBKMB_L4_axons_queries_2017_b'};
% skeletonFolders = cellfun(@(x)[scratchFolder x filesep], skeletonFolders, 'uni', 0);
%%%% old dendrite round
% scratchFolder = '/tmpscratch/kboerg/';
% % These are the downloaded queries (see wK projects: L4_focus_flight-1 & 2, L4_focus_flight-reseed-1 & 2 (2 contains 2nd and 3rd round of reseeding now))
% % And new queries from after switching to new (agglomerate based) query analysis: L4_focus_flight-new-1
% skeletonFolders = {'/tmpscratch/kboerg/MBKMB_L4_dendrite_queries_2017_c_nmls/','/tmpscratch/kboerg/MBKMB_L4_dendrite_queries_2017_c_inverse_nmls/'};
% skeletonFolders = cellfun(@(x)[scratchFolder x filesep], skeletonFolders, 'uni', 0);
%%%%% this round
% skeletonFolders = {'/tmpscratch/scchr/AxonEndings/axonQueryResults/CS_MB_L4_AxonLeftQueries_nmls/'};


% Load segment meta data & graph
load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat');
graph = connectEM.loadAllSegmentationData(p);

if isaxon
    edgesGTall = load('/gaba/scratch/mberning/edgesGTall.mat');
    edgesGTall = edgesGTall.edgesGTall;
    functionH = @connectEM.queryAnalysisSub;
    inputCell = cellfun(@(x){x}, num2cell(1 : 500), 'uni', 0);
    cluster = Cluster.getCluster( ...
        '-pe openmp 1', ...
        '-p 0', ...
        '-l h_vmem=24G', ...
        '-l s_rt=23:50:00', ...
        '-l h_rt=24:00:00');
    job = Cluster.startJob( functionH, inputCell, ...
        'name', 'query', ...
        'cluster', cluster,'sharedInputs',{agglos,edgesGTall,outputFolder});
    Cluster.waitForJob(job);
    for startidx = 1  : 500
        temp = load(fullfile(outputFolder,['output' num2str(startidx)]),'usededges');
        usededges(startidx:500:size(agglos,1)) = temp.usededges(startidx:500:size(agglos,1));
    end
    for idx = 1 : length(agglos)
        assert(isequal(Graph.findConnectedComponents(edgesGTall(usededges{idx},:)),{sort(agglos(idx).nodes(:,4))}));
    end
    usededges2=cell2mat(usededges(~cellfun(@isempty,usededges))');
    usededges3 =ismember(graph.edges,edgesGTall(usededges2,:),'rows');
    
    borderNm = repmat([2000; -2000], 1, 3);
    borderVoxel = round(bsxfun(@times, 1./p.raw.voxelSize, borderNm));
    bboxSmall = p.bbox + borderVoxel';
    borderMeta = load([p.saveFolder 'globalBorder.mat'], 'borderCoM');
    borderMeta.borderCoM(end+1,:)=nan(1,3);
    borderidxs = graph.borderIdx(usededges3);
    borderidxs(isnan(borderidxs))=length(borderMeta.borderCoM);
    
    borderPositions = double(borderMeta.borderCoM(borderidxs,:));
    
    outsideBbox = ~(all(bsxfun(@gt, borderPositions, bboxSmall(:, 1)'), 2) & ...
        all(bsxfun(@lt, borderPositions, bboxSmall(:, 2)'), 2));
    usededges4=usededges3;
    usededges4(graph.prob(usededges3)<0.98&outsideBbox)  = false;
%     usededges4 = cat(1,usededges4,repmat((1:numel(agglos))',1,2)); % add
%     self edges not tested yet
    aggloFiltered = Graph.findConnectedComponents(graph.edges(usededges4,:));
    aggloFiltered = [aggloFiltered; num2cell(setdiff(cell2mat(agglos), cell2mat(aggloFiltered)))];
    calculateLength = @(x)max(pdist(bsxfun(@times, double(borderMeta.borderCoM(x, :)), p.raw.voxelSize)));
    filternan = @(x)x(~isnan(x));
    axonLength = NaN(length(aggloFiltered),1);
    for idx = 1 : length(aggloFiltered)
        axonLength(idx) = max([-1, calculateLength(filternan(cell2mat(graph.neighBorderIdx(aggloFiltered{idx}))))]);
    end
    agglos=aggloFiltered(axonLength>5000);
    
end

%%
%------------------------------------------------------------------------------------------
% load('/tmpscratch/mberning/superagglos.mat')
% for i=1:length(superagglos)
%     agglos{i,1} = superagglos(i).nodes(:,4);
% end
% agglos = cellfun(@(x)x(~isnan(x)),agglos,'uni',0);

% flightPaths_old = structfun(@(x)x.nodes(isnan(x.nodes(:,4)), 1:3),superagglos,'uni',0);
% for i=1:length(superagglos)
%     flightPaths_old{i,1} = superagglos(i).nodes(isnan(superagglos(i).nodes(:,4)), 1:3);
% end



% Lookup segment ids of nodes+neighbours of nmls in all folders defined above
[ff.segIds, ff.neighbours, ff.filenames, ff.nodes, ff.startNode, ff.comments] = connectEM.lookupNmlMulti(p, skeletonFolders, false);
display([num2str(sum(~cellfun(@isempty,ff.comments))) '/' num2str(numel(ff.comments)) ' queries contain comment and will not be used']);
tabulate(cellfun(@(x)x{1}{1}, cellfun(@(x)regexp(x, 'content="(.*)"', 'tokens'), ...
    cat(1, ff.comments{~cellfun(@isempty, ff.comments)}), 'uni', 0), 'uni', 0))
ff = structfun(@(x)x(cellfun(@isempty,ff.comments)), ff, 'uni', 0);
% ~600 queries do not have a start node, not sure why (maybe the ones with more than one tree), maybe check later
ff = structfun(@(x)x(~cellfun(@isempty, ff.startNode)), ff, 'uni', 0);

save(fullfile(outputFolder,'ff_structAggloQueries.mat'), 'ff');
% load('/tmpscratch/scchr/AxonEndings/axonQueryResults/ff_struct_CS_MB_L4_AxonLeftQueries.mat', 'ff');


% Calculate overlap of all queries with segments
[uniqueSegments, neighboursStartNode, ~, ~] = cellfun(@connectEM.queryAnalysis, ...
    ff.segIds, ff.neighbours, ff.nodes', ff.startNode', 'uni', 0);
% Determine all overlaps of agglomerations with given queries
[~, queryOverlap] = connectEM.queryAgglomerationOverlap(agglos, segmentsLeftover, uniqueSegments, neighboursStartNode);
% Make decision(s), here evidence/occurence threshold is applied
% Always one (or none if evidence below 14, 1/2 node) start eqClass
startAgglo = arrayfun(@(x)x.eqClasses(x.occurences > 13), queryOverlap.start, 'uni', 0);
% Exclude all queries that do not have a clear starting point
idxNoClearStart = cellfun('isempty', startAgglo);
% Multiple ends (all above 53vx evidence, corresponds to 2 full nodes)
endAgglo = arrayfun(@(x)x.eqClasses(x.occurences > 53), queryOverlap.ends, 'uni', 0);
% Exclude startAgglo from endAgglo (as we do not want to count self-attachment)
endAgglo = cellfun(@(x,y)setdiff(x,y), endAgglo, startAgglo, 'uni', 0);
% Exclude all queries that do not have (at least one) clear end
idxNoClearEnd = cellfun('isempty', endAgglo);
% 18.5% of queries excluded overall due to missing start or end (or both)
idxGood = ~(idxNoClearStart | idxNoClearEnd);
% Display some statistics
display([num2str(sum(idxNoClearStart)./numel(idxNoClearStart)*100, '%.2f') '% of remaining queries have no clear start']);
display([num2str(sum(idxNoClearEnd)./numel(idxNoClearEnd)*100, '%.2f') '% of remaining queries have no clear end']);
display([num2str(sum(idxGood)./numel(idxGood)*100, '%.2f') '% of remaining queries have clear start and ending']);
display([num2str(numel(cat(2, endAgglo{idxGood}))) ' attachments made by ' num2str(sum(idxGood)) ' queries']);
% Find CC of eqClasses to be joined including single eqClasses with or
% without dangling query
edgesCC = cellfun(@(x,y)combnk([x y], 2), startAgglo(idxGood), endAgglo(idxGood), 'uni', 0);
edgesCC = cat(1,edgesCC{:});
edgesCC = sort(edgesCC, 2);
eqClassCC = Graph.findConnectedComponents(edgesCC, true, true);
sizeEqClassCC = sort(cellfun(@numel, eqClassCC), 'descend');
eqClassCCfull = [eqClassCC; num2cell(setdiff(1 : length(agglos), cell2mat(eqClassCC)))'];
display(sizeEqClassCC(1:10));

results.startAgglo = startAgglo;
results.endAgglo = endAgglo;
results.ff = ff;
results.idxGood = idxGood;

% iterate over super agglos
numstr = '25';

for startingidx = 1:500
    for idx_agglo = startingidx : 500 : length(eqClassCCfull)
        currentEC =eqClassCCfull{idx_agglo};
        mkdir(fullfile(outputFolder,'/visX', [numstr, '_' num2str(floor(idx_agglo/100)) '/']));
        [nodes2, edges2, segIds2] = makeSkelForChiasmataDetectionSub(currentEC, agglos, {results}, false); %!!!!!!!!!!!!!!!!!!!
        assert(length(Graph.findConnectedComponents(edges2))<=1);
        % Save state as used for connectEM.detectChiasmata of nodes, edges and segIds as arrays (table and graph representation will be generated during collection)
        saveFile = fullfile(outputFolder,['visX',numstr ,'_' num2str(floor(idx_agglo/100))], ['visX', numstr, '_' num2str(idx_agglo)], 'superagglo.mat');
        mkdir(fileparts(saveFile));
        save(saveFile, 'nodes2', 'edges2', 'segIds2');
    end
end

connectEM.collectSuperagglos(outputFolder,prefix,outputFolder)

% % Generate new merged classes
% eqClassCCfull = [eqClassCC; num2cell(setdiff(1 : length(agglos), cell2mat(eqClassCC)))'];
% agglosNew = cellfun(@(x){cell2mat(agglos(x))}, eqClassCCfull);
% % Load sorted out (and not used here) small agglomerates and add for Alessandro
% load('/gaba/scratch/mberning/axonQueryGeneration/axonsSmall.mat', 'axonsSmall');
% agglosPostQuery = cat(1, agglosNew, agglosSmall);
% save([outputFolder 'postQueryAgglomerates.mat'], 'agglosPostQuery');

%Save workspace
save([outputFolder 'postQueryState.mat']);
