function nBPs = getBranchPoints(p,superagglos,method,show,useCluster)
if ~exist('show','var') || isempty(show)
    show = 0;
end
if ~exist('useCluster','var') || isempty(useCluster)
    useCluster = 0;
end
if ~exist('method','var') || isempty(method)
    method = 'thin';
end
superagglos = rmfield(superagglos,setdiff(fieldnames(superagglos),{'nodes','edges'})); % remove all fields except nodes/edges

somaAgglos = connectEM.getSomaAgglos(fullfile(p.saveFolder,'aggloState','somas_with_merged_somas.mat'),'all');

sizeInnerBBox = 2000;  % in nm. size of bbox..
innerbbox = [(-sizeInnerBBox/2 ./ p.raw.voxelSize);(sizeInnerBBox/2 ./ p.raw.voxelSize)]';
sizeOuterBBox = 10000;  % in nm. size of bbox..
outerbbox = [(-sizeOuterBBox/2 ./ p.raw.voxelSize);(sizeOuterBBox/2 ./ p.raw.voxelSize)]';
nBPs = zeros(numel(superagglos),1);

aggloWOsoma = cell(numel(superagglos),1);
tic
for s = 1:numel(superagglos)
    thisagglo = Superagglos.removeDuplicates(superagglos(s));
    if useCluster
        % cut out soma
        aggloWOsoma{s} = Superagglos.removeSegIdsFromAgglos(thisagglo,cell2mat(somaAgglos),1);
    else
        nBPs(s) = getBranchPoint(Superagglos.removeSegIdsFromAgglos(thisagglo,cell2mat(somaAgglos),1),method,show,innerbbox,outerbbox);
        Util.progressBar(s,numel(superagglos));
    end
end
cluster = Cluster.getCluster( ...
    '-l h_vmem=12G', ...
    '-l s_rt=3:59:00', ...
    '-l h_rt=4:00:00');

job = Cluster.startJob( ...
    @getBranchPoint, aggloWOsoma, ...
    'sharedInputs', {{method},{show},{innerbbox},{outerbbox}}, ...
    'sharedInputsLocation', 2:5, ...
    'cluster', cluster,'name','branchpoints');
Cluster.waitForJob(job);
nBPs = fetchOutputs(job);
nBPs = cell2mat(nBPs);

function nBPs = getBranchPoint(aggloWOsoma,method,show,innerbbox,outerbbox)
nBPs = 0;

if size(aggloWOsoma.nodes,1) <= 2   % BP can not exist with only 2 nodes or less
    return
end
BPcoords = zeros(0,3);
switch method
    case 'sphere'
        BPs = [];
        % go through all processes left after cutting out soma
        for n = 1:size(aggloWOsoma.nodes,1)
            % constrain superagglos to everything in the outerBbox around the current
            % node and choose that connected agglo which has the current node
            % in it
            nodesToRemove = ~(all(bsxfun(@ge,aggloWOsoma.nodes(:,1:3),aggloWOsoma.nodes(n,1:3)+outerbbox(:,1)'),2) & all(bsxfun(@le,aggloWOsoma.nodes(:,1:3),aggloWOsoma.nodes(n,1:3)+outerbbox(:,2)'),2));
            thisAgglo = Superagglos.removeNodesFromAgglo(aggloWOsoma,nodesToRemove); % remove all nodes outside the outer bbox
            thisAgglo = thisAgglo(arrayfun(@(x) any(ismember(x.nodes(:,1:3),aggloWOsoma.nodes(n,1:3),'rows')),thisAgglo)); % get the superagglos which belongs to the current node
            % cut out innerBbox around the current node
            nodesToRemove = (all(bsxfun(@ge,thisAgglo.nodes(:,1:3),aggloWOsoma.nodes(n,1:3)+innerbbox(:,1)'),2) & all(bsxfun(@le,thisAgglo.nodes(:,1:3),aggloWOsoma.nodes(n,1:3)+innerbbox(:,2)'),2));
            splitAgglo = Superagglos.removeNodesFromAgglo(thisAgglo,nodesToRemove); % remove all nodes that are not connected anymore
            % if there are more than 2 left connected components which are
            % longer than 500 nm, remember node as possible branchpoint
            if sum(Superagglos.calculateTotalPathLength(splitAgglo,[11.24,11.24,28]) > 500) > 2
                BPs = cat(1,BPs,n);
            end
        end
        % cluster all found BP candidates, as there are multiple at each
        % location
        BPClusters = Graph.findConnectedComponents(aggloWOsoma.edges(any(ismember(aggloWOsoma.edges,BPs),2),:),1);
        % take one coordinate per cluster and save it. save also number of
        % clusters ( = BPs)
        BPcoords = cat(1,BPcoords,aggloWOsoma.nodes(cellfun(@(x) x(1),BPClusters),1:3));
        nBPs = nBPs + numel(BPClusters);
        
        BPcoords = bsxfun(@times,BPcoords,[11.24,11.24,28]);
        if show
            fprintf('Cutting out approach found %d BPs in superagglos %d.\n',nBPs,s);
        end
    case 'thin'
        %% other approach: thin out skeleton
        initialNumCCs = numel(Graph.findConnectedComponents(cat(1,aggloWOsoma.edges,repmat(unique(aggloWOsoma.edges(:)),1,2)),0,1));
        % sort edges to have the highly connected ones checked first
        [uEdges,~,uIdx] = unique(aggloWOsoma.edges);
        counts = histc(aggloWOsoma.edges(:),uEdges);
        if ~any(counts>2) % linear agglo, skip
            return
        end
        [~,idx] = sort(sum(reshape(counts(uIdx),[],2),2));
        aggloWOsoma.edges = aggloWOsoma.edges(idx,:);
        distThr = 1000; % nm
        vipEdgeCounter = 0;
        while 1
            % remove last edge from edge list and check if aggloWOsoma does not fall
            % apart
            edges = aggloWOsoma.edges(1:end-1,:);
            equivalenceClass = Graph.findConnectedComponents(cat(1,edges,repmat(unique(aggloWOsoma.edges(:)),1,2)),0,1);
            if size(edges,1) <= vipEdgeCounter
                break
            end
            % check if agglo did not split or if it split but the only one
            % segment was splitted apart which was not more distant than 1 um.
            if numel(equivalenceClass) == initialNumCCs || ( numel(equivalenceClass) == initialNumCCs+1 && any(cellfun(@numel,equivalenceClass)==1) && sqrt(sum((diff(aggloWOsoma.nodes(aggloWOsoma.edges(end,:),1:3),1,1).*[11.24 11.24 28]).^2)) < distThr)
                aggloWOsoma.edges = edges;  % overwrite the edges with the new edge list (missing the unneeded edge)
            else
                % if edges is needed, move it to front and remember number of
                % important edges
                aggloWOsoma.edges = circshift(aggloWOsoma.edges,1,1);
                vipEdgeCounter = vipEdgeCounter + 1;
            end
        end
        aggloWOsoma = Superagglos.removeNodesFromAgglo(aggloWOsoma,setdiff(1:size(aggloWOsoma.nodes,1),aggloWOsoma.edges(:)),1); % remove all nodes that are not connected anymore
        if isempty(aggloWOsoma)   % BP can not exist with only 2 nodes or less
            return
        end
        aggloWOsoma = Superagglos.removeNodesFromAgglo(aggloWOsoma,histc(aggloWOsoma.edges(:),1:size(aggloWOsoma.nodes,1)) == 1); % remove all endings. This of course also shortens the aggloWOsoma bit but should remove most of the untrue endings
        BPcoords = zeros(0,3);
        for a = 1:numel(aggloWOsoma)
            BPcands = find(histc(aggloWOsoma(a).edges(:),1:size(aggloWOsoma(a).nodes,1)) > 2);
            BPs = BPcands(arrayfun(@(x) sum(Superagglos.calculateTotalPathLength(Superagglos.removeNodesFromAgglo(aggloWOsoma(a),x),[11.24,11.24,28]) > 500)> 2 ,BPcands));
            
            %     aggloWOsoma = Superagglos.removeSegIdsFromAgglos(aggloWOsoma,cell2mat(somaAgglos));
            %     BPcoords = cell2mat(arrayfun(@(x) x.nodes(histc(x.edges(:),1:size(x.nodes,1)) > 2,1:3),aggloWOsoma,'uni',0)); % find putative branch points
            
            BPcoords = cat(1,BPcoords,bsxfun(@times,aggloWOsoma(a).nodes(BPs,1:3),[11.24,11.24,28]));
            nBPs = nBPs+numel(BPs);
        end
        if show
            fprintf('Thinning out approach found %d BPs in superagglo %d.\n',nBPs,s);
        end
    otherwise
        error('Method unknown');
end

if show
    figure;hold all
    skel = Superagglos.toSkel(aggloWOsoma);
    %         skel = Superagglos.toSkel(aggloWOsoma);
    skel.plot;
    axis equal
    %         scatter3(aggloWOsoma.nodes(aggloWOsoma.edges(end,:),1)*11.24,aggloWOsoma.nodes(aggloWOsoma.edges(end,:),2)*11.24,aggloWOsoma.nodes(aggloWOsoma.edges(end,:),3)*28,150,'r','o','filled')
    scatter3(BPcoords(:,1),BPcoords(:,2),BPcoords(:,3),150,'r','o','filled')
end


