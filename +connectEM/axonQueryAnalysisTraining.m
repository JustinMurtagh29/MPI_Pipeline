addpath(genpath('/gaba/u/kboerg/code/auxiliaryMethods/'));
scratchFolder = '/tmpscratch/kboerg/'
skeletonFolders = {'HW_L4_GroundTruthForFlightTraining_nmls'};
skeletonFolders = cellfun(@(x)[scratchFolder x filesep], skeletonFolders, 'uni', 0);
% Lookup segment ids of nodes+neighbours of nmls in all folders defined above
load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat');
temp = load('/gaba/scratch/mberning/axonQueryGeneration/beforeQueryGeneration.mat', 'axonsNew');
axons = temp.axonsNew;
clear temp
disp('calculating ground truth...')

[ff.segIds, ff.neighbours, ff.filenames, ff.nodes, ff.startNode, ff.comments] = connectEM.lookupNmlMulti(p, skeletonFolders, false);
ff = structfun(@(x)x(cellfun(@isempty,ff.comments)), ff, 'uni', 0);
ff = structfun(@(x)x(~cellfun(@isempty, ff.startNode)), ff, 'uni', 0);

% Calculate overlap of all queries with segments
[uniqueSegments, neighboursStartNode, nodesExcludedIdx, startNodeIdx] = cellfun(@connectEM.queryAnalysis, ...
        ff.segIds, ff.neighbours, ff.nodes', ff.startNode', 'uni', 0);
% Determine all overlaps of agglomerations with given queries
[partition, queryOverlap] = connectEM.queryAgglomerationOverlap(axons, [], uniqueSegments, neighboursStartNode);
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

% Look at redundant control tasks
liste = ff.filenames;
last = @(x)x{end};
fifth = @(x)x{5};
sixth = @(x)x{6};
getTask = @(x){fifth(strsplit(last(strsplit(x,'/')),'_'))};
getUser = @(x){sixth(strsplit(last(strsplit(x,'/')),'_'))};

taskStrings = unique(cellfun(getTask, liste))';
taskStrings{1,2} = [];
for idx = 1 : length(liste)
    %if mod(idx, 100) == 0
    %    idx
    %end
    rowidx = strcmp(taskStrings(:, 1), getTask(liste{idx}));
    taskStrings{rowidx, 2} = [taskStrings{rowidx, 2}, idx];
end
disp('ground truth calculation finished.')

scratchFolder = '/u/hwissler/temp/'
skeletonFolders = {'HW_flight_training_2'}
skeletonFolders = cellfun(@(x)[scratchFolder x filesep], skeletonFolders, 'uni', 0);

disp('preparing cells in given folder and compare with ground truth...')
% Lookup segment ids of nodes+neighbours of nmls in all folders defined above
[ff2.segIds, ff2.neighbours, ff2.filenames, ff2.nodes, ff2.startNode, ff2.comments] = connectEM.lookupNmlMulti(p, skeletonFolders, false);
%ff2 = structfun(@(x)x(cellfun(@isempty,ff2.comments)), ff2, 'uni', 0); %löscht files mit comments
ff2 = structfun(@(x)x(~cellfun(@isempty, ff2.startNode)), ff2, 'uni', 0); %löscht files mit fehlender startnode

[uniqueSegments2, neighboursStartNode2, nodesExcludedIdx2, startNodeIdx2] = cellfun(@connectEM.queryAnalysis, ...
        ff2.segIds, ff2.neighbours, ff2.nodes', ff2.startNode', 'uni', 0);

% Determine all overlaps of agglomerations with given queries
[partition2, queryOverlap2] = connectEM.queryAgglomerationOverlap(axons, [], uniqueSegments2, neighboursStartNode2);
% Make decision(s), here evidence/occurence threshold is applied
% Always one (or none if evidence below 14, 1/2 node) start eqClass
startAgglo2 = arrayfun(@(x)x.eqClasses(x.occurences > 13), queryOverlap2.start, 'uni', 0);
% Exclude all queries that do not have a clear starting point
idxNoClearStart2 = cellfun('isempty', startAgglo2);
% Multiple ends (all above 53vx evidence, corresponds to 2 full nodes)
endAgglo2 = arrayfun(@(x)x.eqClasses(x.occurences > 53), queryOverlap2.ends, 'uni', 0);
% Exclude startAgglo from endAgglo (as we do not want to count self-attachment)
endAgglo2 = cellfun(@(x,y)setdiff(x,y), endAgglo2, startAgglo2, 'uni', 0);
% Exclude all queries that do not have (at least one) clear end
idxNoClearEnd2 = cellfun('isempty', endAgglo2);
% 18.5% of queries excluded overall due to missing start or end (or both)
idxGood2 = ~(idxNoClearStart2 | idxNoClearEnd2);

disp('calculating results and prepare to show in table')
unames = unique(cellfun(getUser,ff2.filenames))';
unames{1, 2} = [];
for idx2 = 1 : length(ff2.startNode)
    idxs = find(ismember(cell2mat(ff.startNode'),ff2.startNode{idx2},'rows'));
    if isempty(idxs)
        %warning('unconnected tracing')
        unames{strcmp(unames(:,1),getUser(ff2.filenames{idx2})),2}(end+1) = -2;
    else
        idxsrow = find(cellfun(@(x)ismember(idxs(1),x),taskStrings(:,2)));
        assert(isequal(taskStrings{idxsrow,2},idxs'));
        tab2start=sortrows(tabulate(cellfun(@mat2str,startAgglo(idxs),'uni',0)),2);
        tab2end=sortrows(tabulate(cellfun(@mat2str,endAgglo(idxs),'uni',0)),2);
        if strcmp(tab2end(end,1),mat2str(endAgglo2{idx2})) && strcmp(tab2start(end,1),mat2str(startAgglo2{idx2}))
            unames{strcmp(unames(:,1),getUser(ff2.filenames{idx2})),2}(end+1) = 1;
        else
            unames{strcmp(unames(:,1),getUser(ff2.filenames{idx2})),2}(end+1) = -1;
        end
    end
end
unames(:,3)=cellfun(@sum,unames(:,2),'uni',0);
unames(:,4)=cellfun(@length,unames(:,2),'uni',0);
unames(:,5)=cellfun(@(x)sum(x==-2),unames(:,2),'uni',0);
unames(:,6)=cellfun(@(x)sum(x==-1),unames(:,2),'uni',0);
unames(:,7)=cellfun(@(x)sum(x==1),unames(:,2),'uni',0);
disp('unames:')
disp('......................allResultsTable, score,   #all, #corrupt, #wrong, #correct, Percentage correct')
unames(:,8)=cellfun(@(x)sum(x==1)/length(x)*100,unames(:,2),'uni',0)
