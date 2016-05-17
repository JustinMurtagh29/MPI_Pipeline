function boutonAgglomerationFull(p)

% Set rng to state of starting matlab
rng default;

%% Load some stuff
% Load parameter from newest pipeline run
load([p.saveFolder 'allParameter.mat']);
% Load global graph representation
graph = load([p.saveFolder 'graph.mat']);
% Load bouton seeds (presynaptic segments as detected by Benedikt)
boutons = load([p.saveFolder 'Boutons.mat']);
% Load center of mass of all segments
load([p.saveFolder 'globalCoMList.mat']);
% Load ground truth already agglomerated segments 
groundTruth = load([p.saveFolder 'agglomeration/20160202T222401.mat']);
% Fix excludeIds to represent @94% 
excludeIds = cat(1, groundTruth.ids{:,7}, groundTruth.excludeIds);

%% Seed generation
% Restrict seeds to the ones we are interested in (within segmented region and not yet agglomerated)
seeds.all = unique(intersect(boutons.boutonIDs,unique(graph.edges(:))));
seeds.notTraced = setdiff(seeds.all, excludeIds);
% Drop those with tile distance 2 or less or less to border (e.g. approx 10 microns on each sice)
bbox = p.bbox + [1024 -1024; 1024 -1024; 512 -512];
idx = all(bsxfun(@gt, boutons.boutonCoMs, bbox(:,1)'),2) & all(bsxfun(@lt, boutons.boutonCoMs, bbox(:,2)'),2);
seeds.restricted = boutons.boutonIDs(idx);
clear idx bbox;

%% Bouton agglomeration
% Agglomerate 1000 random seeds chosen above @90% and @40 neighbours thresholds
cluster = createLocalParcluster(); 
parpool(cluster,8);
parfor i=1:length(seeds.restricted)
    [collectedIds{i}, probabilities{i}, excluded(i)] = agglomerateSG6(graph, com, seeds.restricted(i), excludeIds, 0.90, 50);
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
tic;
for i=1:length(seeds.restricted)
    idxInVec = find(ccVector == seeds.restricted(i));
    chosenSeedIdxInCC(i) = find(ccLength >= idxInVec, 1, 'first');
    Util.progressBar(i, length(seeds.restricted));
end
% Decide which positions (or edges for games) to querry
tic;
for i=1:length(cc)
    [q.pos{i}, q.dir{i}, q.edges{i}, q.varianceExplained{i}] = determineQuerryLocation(graph, com, cc{i});
    Util.progressBar(i, length(cc))
end
clear i edgesAgglo ccVector ccLength idxInVec skelToWrite theseSeedIdx skelName;

%% Sort out agglomerations based on exclusion list
nrSeeds = cellfun(@(x)length(intersect(seeds.notTraced,x)), cc);
excludedCC = false(length(cc),1);
tic;
for i=1:length(cc)
    theseSeedIdx = chosenSeedIdxInCC == i;
    reasons = cat(2,excluded(theseSeedIdx).reasons);
    nrExcludedThisCC = sum(cellfun(@(x)strcmp(x, 'exclusion list') ~= 0, reasons));
    if nrExcludedThisCC > 1
        excludedCC(i) = true;
    end
    Util.progressBar(i, length(cc));
end

% Add agglomerations were PC1 explains less than 70%
excludedCC = excludedCC | (cellfun(@(x)x(1), q.varianceExplained)' < 0.7 & cellfun(@length, cc) > 40);
for i=1:length(cc)
    theseSeedIdx = chosenSeedIdxInCC == i;
    seedIds = intersect(seeds.notTraced, cc{i});
    skelToWrite = writeSkeletonFromAggloWBP(graph, com, cc(i), excluded(theseSeedIdx), seedIds, q.edges{i}(:,1), ['CC' num2str(i, '%.3i')]);
    if excludedCC(i)
        skelName = [p.saveFolder 'axonsSorted/bad' num2str(i, '%.3i') '.nml'];
    else 
        skelName = [p.saveFolder 'axonsSorted/good' num2str(i, '%.3i') '.nml'];
    end
    evalc('writeNml(skelName, skelToWrite)');
end

% 'Write' problems to flight mode webKNOSSOS
extend = round(1000 ./ [11.24 11.24 28]);
fid = fopen([p.saveFolder 'webKnossos.txt'], 'w');
for i=1:length(q.pos)
    if ~excludedCC(i)
        for j=1:size(q.pos{i},1)
            [phi, theta, psi] = calculateEulerAngles(-q.dir{i}(j,:));
            minPos = q.pos{i}(j,:) - extend; 
            maxPos = q.pos{i}(j,:) + extend;
            linkString = ['2012-09-28_ex145_07x2,focus_flight,focus_flight,1,' ...
                num2str(q.pos{i}(j,1)) ',' num2str(q.pos{i}(j,2)) ',' num2str(q.pos{i}(j,3)) ',' ...
                num2str(phi) ',' num2str(theta) ',' num2str(psi) ',100,1,Tracing crew,' ...
                num2str(minPos(1)) ',' num2str(minPos(2)) ',' num2str(minPos(3)) ',' ...
                num2str(maxPos(1)) ',' num2str(maxPos(2)) ',' num2str(maxPos(3)) ',' 'L4_focus_flight-1'];
            fprintf(fid, '%s\n', linkString);
            q.angles{i}(j,:) = [phi theta psi];
        end
    end
end
fclose(fid);

% Save variables 
variables = who;
% Exclude large ones loaded in beginning
variables(strcmp(variables, 'graph')) = [];
variables(strcmp(variables, 'com')) = [];
variables(strcmp(variables, 'groundTruth')) = [];
save([p.saveFolder ' datestr(clock, 30) '_boutonAggloFull.mat'], variables{:});

end
