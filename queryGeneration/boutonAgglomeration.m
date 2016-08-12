
function boutonAgglomeration(p)
% Set rng to state of starting matlab
rng default;
% First test with global agglomerations (set some parameter here)
n_axon_agglo = 50;
p_axon_agglo = .95;
p_nuclei_agglo = 7;

%% Load some stuff
% Load global graph representation
graphStruct = load([p.saveFolder 'graph.mat']);
% Load bouton seeds (presynaptic segments as detected by Benedikt)
m= load([p.saveFolder 'Boutons.mat']);
boutons=m.boutons;
% Load center of mass of all segments
load([p.saveFolder 'globalCoMList.mat']);
com=globalCoMList;
% Load inital agglomeration of vessel, nuclei border
m=load([p.saveFolder 'nucleiSegIds.mat']);
nucleiSegIds = m.segIds;
% Output folder 
outputFolder = ['/gaba/u/ganja/results_P14_L4_Corrected/Skeletons' datestr(clock,30) '/'];
mkdir(outputFolder);
mkdir([outputFolder 'axons']);
mkdir([outputFolder 'axonsSorted']);
mkdir([outputFolder 'movies']);
mkdir([outputFolder 'nucleiAgglo']);

%% Nuclei and Vessel Agglomeration cleanup
% Clean up agglomeration (remove earlier agglomerations & split up vessel according to cc)
%agglo = rmfield(agglo, {'borderGrown' 'nucleiGrown'});
% Find overlaps between components (exclude borders from check, as those CAN be duplicates of either) 
%[betweenLabels, withinLabels] = checkAggloConsistency(rmfield(agglo, 'borderMerged'));
% Remove detected overlaps between nuclei and vessel detection & within nuclei & check again
%fH = @(x)setdiff(x,cat(1,betweenLabels{:}));
%aggloNew.nuclei = cellfun(fH, agglo.nuclei, 'UniformOutput', false);
%aggloNew.vessel = cellfun(fH, agglo.vessel, 'UniformOutput', false);
%fH = @(x)setdiff(x,cat(1,withinLabels{:}));
%aggloNew.nuclei = cellfun(fH, aggloNew.nuclei, 'UniformOutput', false);
%[betweenLabels, withinLabels] = checkAggloConsistency(aggloNew);
% Checked, after this no more overlaps, now remove empty (<= 5 nodes) components and clean up
%fN = fieldnames(aggloNew);
%for i=1:length(fN)
%    for j=1:length(aggloNew.(fN{i}))
%        temp = findCCaccordingToGraph(graphStruct, aggloNew.(fN{i}){j});
%        sizeCC = cellfun('length', temp);
%        idx = sizeCC > 5;
%        aggloNew.(fN{i}){j} = temp(idx);
%    end
%end
%aggloNew.nuclei = cat(1,aggloNew.nuclei{:});
%aggloNew.vessel = cat(1,aggloNew.vessel{:});
%aggloNew.border = agglo.borderMerged;
%clear fH fN agglo betweenLabels withinLabels i j idx ans temp sizeCC pT;
% Write vessel and nuclei to nml for first checks
%writeNml([outputFolder 'vessel.nml'], writeSkeletonFromAgglo(graphStruct, com, aggloNew.vessel, 'Vessel '));
%writeNml([outputFolder 'nuclei.nml'], writeSkeletonFromAgglo(graphStruct, com, aggloNew.nuclei, 'Nuclei '));

%% Now agglomerate detected nuclei in 1% decrement down to 94%
% Detect overlap between agglomerations at each step and exclude
prob = [1.00:-0.01:.94];
%excludeIds = cat(1, aggloNew.nuclei{:}, aggloNew.vessel{:});

aggloNew.nuclei = nucleiSegIds;
excludeIds =  cat(1, aggloNew.nuclei{:});
for p_i=1:length(prob)
    for i=1:length(aggloNew.nuclei)
        if p_i == 1
            ids{i,p_i} = aggloNew.nuclei{i};
            nrAgglo(i,p_i) = 0;
        else
            [ids{i,p_i}, nrAgglo(i,p_i)] = agglomerateSG_simple(graphStruct, ids{i, p_i-1}, prob(p_i), excludeIds);
        end
    end
    % Check for & remove overlaps if any are found
    if p_i > 1
        newIds = cellfun(@setdiff, ids(:,p_i), ids(:,p_i-1), 'UniformOutput', false);
        idx = cellfun('isempty', newIds);
        newIds(idx) = [];
        newExcludeIds{p_i} = unique(cat(1,newIds{:}));
        excludeIds = cat(1,excludeIds,newExcludeIds{p_i});
        overlap = cell(length(newIds));
        for k=1:length(newIds)
            for j=k+1:length(newIds)
                temp = intersect(newIds{k},newIds{j});
                if ~isempty(temp)
                    overlap{k,j} = temp;
                end
            end
        end
        overlapIds{p_i} = cat(1, overlap{:});
        keepIdx = cellfun(@(x)~ismember(x, overlapIds{p_i}), ids(:,p_i),'UniformOutput', false);
        ids(:,p_i) = cellfun(@(x,y)x(y), ids(:,p_i), keepIdx, 'UniformOutput', false);
        nrExcluded(p_i) = sum(cellfun(@(x)sum(~x),keepIdx));
    end
end
clear i p_i j k overlap keepIdx newIds idx temp;
save([outputFolder 'nucleiAgglomeration.mat'], 'ids', 'nrAgglo', 'overlapIds', 'nrExcluded');

%% Visualize results of nuclei agglomeration
% Look at the agglomerations (minimal spanning tree on COM distances of agglomerated segments)
%for i=1:size(ids,1)
%    skel = writeSkeletonFromAggloCombinedMerged(graphStruct, com, ids(i,:), 'Probability: ', prob);
%    skelName = [outputFolder 'nucleiAgglo/nucleiAgglo' num2str(i, '%.4i') '.nml'];
%    evalc('writeNml(skelName, skel)');
%end
% Looking at agglomerations in these skeletons and some statistics on file saved above, i will now take p_i = 7
% Look at excluded ids from nuclei agglomeration in wK for up to some probability threshold
%leftoverID = cat(1,overlapIds{1:p_nuclei_agglo});
%leftoverCC = findCCaccordingToGraph(graph, leftoverID);
%skel = writeSkeletonFromAgglo(graph, com, leftoverCC);
%writeNml([outputFolder 'leftovers.nml'], skel);
% This time final skeleton @94% to use for exclusion list
%writeNml([outputFolder 'nucleiGrown.nml'], writeSkeletonFromAgglo(graph, com, ids(:,p_nuclei_agglo), 'Nuclei '));
%clear skel leftoverID leftoverCC i;

%% Seed generation
% Fix excludeIds to represent @94% 
%excludeIds = cat(1, ids{:,p_nuclei_agglo}, aggloNew.vessel{:}, overlapIds{1:p_nuclei_agglo});
excludeIds = cat(1, ids{:,p_nuclei_agglo}, overlapIds{1:p_nuclei_agglo});
% Restrict seeds to the ones we are interested in (within segmented region and not yet agglomerated)
seeds.all = unique(intersect(boutons.boutonIDs,unique(graphStruct.edges(:))));
seeds.notTraced = setdiff(seeds.all, excludeIds);
% Drop those with tile distance 2 or less or less to border (e.g. approx 10 microns on each sice)
bbox = p.bbox + [512 -512;512 -512;256 -256];
idx = all(bsxfun(@gt, boutons.boutonCoMs, bbox(:,1)'),2) & all(bsxfun(@lt, boutons.boutonCoMs, bbox(:,2)'),2);
seeds.restricted = boutons.boutonIDs(idx);
% Choose 1000 randomly
seeds.chosen = seeds.restricted;
clear idx bbox;

%% Bouton agglomeration
% Agglomerate 1000 random seeds chosen above @80% and @40 neighbours thresholds
for i=1:length(seeds.chosen)
    [collectedIds{i}, probabilities{i}, excluded(i)] = agglomerateSG6(graphStruct, com, seeds.chosen(i), excludeIds, p_axon_agglo, n_axon_agglo);
end
% Construct edge list for finding CC of agglomerations 
edgesAgglo = cellfun(@(x)combnk(x,2), collectedIds, 'UniformOutput', false);
edgesAgglo = cat(1, edgesAgglo{:});
edgesAgglo = sort(edgesAgglo,2);
edgesAgglo = unique(edgesAgglo, 'rows');
% Find connected components of agglomerated ids
cc = Graph.findConnectedComponents(edgesAgglo, false, true);
% Find connected component that each CHOSEN seed belongs to
ccVector = cat(1, cc{:});
ccLength = cumsum(cellfun('length', cc)); % apparently string functions are faster, not important here
for i=1:length(seeds.chosen)
    idxInVec = find(ccVector == seeds.chosen(i));
    chosenSeedIdxInCC(i) = find(ccLength >= idxInVec, 1, 'first');
end
% Decide which positions (or edges for games) to querry
for i=1:length(cc)
    [q.pos{i}, q.dir{i}, q.edges{i}, q.varianceExplained{i}] = determineQuerryLocation(graphStruct, com, cc{i});
end


% Visualize all chosen seed agglomerations as skeletons in webKNOSSOS
for i=1:length(cc)
    theseSeedIdx = chosenSeedIdxInCC == i;
    seedIds = intersect(seeds.notTraced, cc{i});
    skelToWrite = writeSkeletonFromAggloWBP(graphStruct, com, cc(i), excluded(theseSeedIdx), seedIds, q.edges{i}(:,1), ['CC' num2str(i, '%.3i')]);
    skelName = [outputFolder 'axons/axons' num2str(i, '%.3i') '.nml'];
    evalc('writeNml(skelName, skelToWrite)');
end

clear i edgesAgglo ccVector ccLength idxInVec skelToWrite theseSeedIdx skelName;

%% Sort out agglomerations based on exclusion list
nrSeeds = cellfun(@(x)length(intersect(seeds.notTraced,x)), cc);
excludedCC = false(length(cc),1);
for i=1:length(cc)
    theseSeedIdx = chosenSeedIdxInCC == i;
    reasons = cat(2,excluded(theseSeedIdx).reasons);
    nrExcludedThisCC = sum(cellfun(@(x)strcmp(x, 'exclusion list') ~= 0, reasons));
    if nrExcludedThisCC > 1
        excludedCC(i) = true;
    end
end
% Add agglomerations were PC1 explains less than 70%
excludedCC = excludedCC | (cellfun(@(x)x(1), q.varianceExplained)' < 0.7 & cellfun(@length, cc) > 40);
for i=1:length(cc)
    theseSeedIdx = chosenSeedIdxInCC == i;
    seedIds = intersect(seeds.notTraced, cc{i});
    skelToWrite = writeSkeletonFromAggloWBP(graphStruct, com, cc(i), excluded(theseSeedIdx), seedIds, q.edges{i}(:,1), ['CC' num2str(i, '%.3i')]);
    if excludedCC(i)
        skelName = [outputFolder 'axonsSorted/bad' num2str(i, '%.3i') '.nml'];
    else 
        skelName = [outputFolder 'axonsSorted/good' num2str(i, '%.3i') '.nml'];
    end
    evalc('writeNml(skelName, skelToWrite)');
end

%{
%% Find 'right' query location, e.g. end of segment along query direction and x,y position
tic;
tempIdx = ismember(graphStruct.edges, sort(cat(1,q.edges{:}),2), 'rows');
tempEdges = graphStruct.edges(tempIdx,:);
tempProb = graphStruct.prob(tempIdx,:);
tempBorderCom = graphStruct.borderCentroid(tempIdx,:);
for i=1:length(q.edges)
    if ~excludedCC(i)
        for j=1:size(q.pos{i},1)
            idx = find(ismember(tempEdges, sort(q.edges{i}(j,:),2), 'rows'));
            if length(idx) > 1
                % In case more than one border, use more likely one
                [~, maxIdx] = max(tempProb(idx));
                idx = idx(maxIdx);
            end
            q.posFixed{i}(j,:) = round(tempBorderCom(idx,:));
        end
    end
end
toc;
clear idx temp* i j;
%}
%{
%% Visualization in movies and generate links for wK
% Visualize querried edges (analog to B4B game 1)
%for i=1:10%length(q.pos)
%    if ~excludedCC(i)
%        for j=1:size(q.pos{i},1)
%            outputFile = [outputFolder 'movies/query_axon' num2str(i, '%.2i') '_pos' num2str(j, '%.2i') '.avi'];
%            arbitraryResliceAgglo(p, graphStruct, com, cc{i}, q.pos{i}(j,:), q.dir{i}(j,:), outputFile);
%        end
%    end
%end
%}
% 'Write' problems to flight mode webKNOSSOS
extend = round(1000 ./ p.raw.voxelSize);
fid = fopen([outputFolder 'webKnossos.txt'], 'w');
for i=1:length(q.pos)
    if ~excludedCC(i)
        for j=1:size(q.pos{i},1)
            [phi, theta, psi] = calculateEulerAngles(-q.dir{i}(j,:));
            minPos = q.pos{i}(j,:) - extend; 
            sizeBBox = 2*extend;
            linkString = ['AG_P14_st003_15-09-2015,focus_flight,focus_flight,1,' ...
                num2str(q.pos{i}(j,1)) ',' num2str(q.pos{i}(j,2)) ',' num2str(q.pos{i}(j,3)) ',' ...
                num2str(phi) ',' num2str(theta) ',' num2str(psi) ',100,1,Tracing crew,' ...
                num2str(minPos(1)) ',' num2str(minPos(2)) ',' num2str(minPos(3)) ',' ...
                num2str(sizeBBox(1)) ',' num2str(sizeBBox(2)) ',' num2str(sizeBBox(3)) ',' 'P14_L4_focus_flight_1'];
            fprintf(fid, '%s\n', linkString);
            q.angles{i}(j,:) = [phi theta psi];
        end
    end
end
fclose(fid);
%{
fid = fopen([outputFolder 'webKnossosFixed.txt'], 'w');
for i=1:length(q.pos)
    if ~excludedCC(i)
        for j=1:size(q.pos{i},1)
            [phi, theta, psi] = calculateEulerAngles(-q.dir{i}(j,:));
            referenceString = ['CC' num2str(i, '%.4i') '_Query' num2str(j, '%.1i')];
            linkString = ['https://webknossos.brain.mpg.de/annotations/Explorational/56bf18a01500005664667830#' ...
                num2str(q.posFixed{i}(j,1)) ',' num2str(q.posFixed{i}(j,2)) ',' num2str(q.posFixed{i}(j,3)) ',1,1.00,' ...
                num2str(phi) ',' num2str(theta) ',' num2str(psi)];
            fprintf(fid, '%s %s\n', referenceString, linkString);
            q.angles{i}(j,:) = [phi theta psi];
        end
    end
end
fclose(fid);

%% Comment section of lablog post
ccFinal = cc(~excludedCC);
pathLength = getPathLengthFromAgglomeration(ccFinal, com);
display(['Parameter: ' num2str(p_axon_agglo, '%3.2f') ' - ' num2str(n_axon_agglo)]);
display(['Average path length: ' num2str(mean(pathLength))]);

% Save variables 
variables = who;
% Exclude large ones loaded in beginning
variables(strcmp(variables, 'graphStruct')) = [];
variables(strcmp(variables, 'com')) = [];
save(['/gaba/u/mberning/' datestr(clock, 30) '_boutonAgglo.mat'], variables{:});
%}
end
