function splitChiasmaSub2(graph,idx_start)
    load('/tmpscratch/kboerg/chiasmarunAugust/preSplitRun');
    p.voxelSize = [11.24 11.24 28];
    for idx = idx_start:50:length(idx_pre)
        superagglo_definition = load(['/tmpscratch/kboerg/visX23/visX19_' num2str(floor(idx_pre(idx)/100)) '/visX19_' num2str(idx_pre(idx)) '/result.mat']);
        %copyfile(['/tmpscratch/kboerg/visX23/visX19_' num2str(floor(idx_pre(idx)/100)) '/visX19_' num2str(idx_pre(idx)) '/skelDetected.nml'],[outputFolder 'skelDetectedOld' num2str(idx) '.nml']);
        % get skels' ids
        skelsegids = connectEM.lookupSkelGT(p, struct('nodes',{{superagglo_definition.output.nodes}}));
        nodesScaled = bsxfun(@times,superagglo_definition.output.nodes,p.voxelSize);
        selected = lookup_tasks(lookup_tasks(:,1)==idx_pre(idx)&idx_pre_pre',2);
        thisEdgesNew = [];
        nodesToAdd = [];
        edgesToDelete=[];
        thisEdgesCol = superagglo_definition.output.edges;
        for idx2 = 1 : length(selected);
            p.sphereRadiusOuter = 10000; % in nm
            p.sphereRadiusInner = 1000; % in nm
            [thisNodes, thisEdges, thisProb] = ...
                connectEM.detectChiasmataPruneToSphere( ...
                nodesScaled, superagglo_definition.output.edges, ...
                ones(size(superagglo_definition.output.edges,1),1), p, superagglo_definition.output.ccCenterIdx(selected(idx2)));
            C = Graph.findConnectedComponents(thisEdges);
            assert(length(C) > 3 && sum(cellfun(@(idx)max(pdist2(thisNodes(idx, :), nodesScaled(superagglo_definition.output.ccCenterIdx(selected(idx2)),:))) > 3000, C))>3);
            nrExits = length(C);
            ff_idx = find(cellfun(@length,strfind(ff.filenames,M(ismember(lookup_tasks,[idx_pre(idx),selected(idx2),1],'rows')))));
            %copyfile(ff.filenames{ff_idx},[outputFolder 'query' num2str(idx) '_' num2str(idx2) '.nml']);
            if nrExits ~= 4
                continue
            end
            getIdsFromCC=@(x)skelsegids{1}(ismember(nodesScaled,thisNodes(x,:),'rows'));
            nonull=@(x)x(x~=0);
            whichones = find(cellfun(@(x)any(ismember(getIdsFromCC(x),nonull(ff.segIds{ff_idx}))),C));
            if length(whichones) ~=2
                continue
            end
            %somehow we lost the node order for the query, here reconstructed with MSP
            Tree = graphminspantree(sparse(squareform(pdist(ff.nodes{ff_idx}))));
            [X,Y]=find(Tree);
            thisEdgesNew = [thisEdgesNew; [X,Y]+size(nodesToAdd,1)+size(superagglo_definition.output.nodes,1)];
            
            for nodeidx = 1 : length(ff.nodes{ff_idx})
                matches = find(ismember(skelsegids{1},ff.segIds{ff_idx}(nodeidx)));
                thisEdgesNew = [thisEdgesNew; matches, repmat(nodeidx+size(nodesToAdd,1)+size(superagglo_definition.output.nodes,1),size(matches))];
            end
            nodesToAdd = [nodesToAdd; ff.nodes{ff_idx}];
            p.sphereRadiusOuter = Inf; % in nm
            p.sphereRadiusInner = 1000; % in nm
            p.voxelSize = [11.24 11.24 28];
            [thisNodes2, thisEdges2, thisProb2] = ...
                connectEM.detectChiasmataPruneToSphere( ...
                nodesScaled, superagglo_definition.output.edges, ...
                ones(size(superagglo_definition.output.edges,1),1), p, superagglo_definition.output.ccCenterIdx(selected(idx2)));
            [~,lookupnodes]=ismember(thisNodes,nodesScaled,'rows');
            [~,lookupnodes2]=ismember(thisNodes2,nodesScaled,'rows');
            edgesNotOuter = sort(setdiff(superagglo_definition.output.edges,setdiff(lookupnodes2(thisEdges2),lookupnodes(thisEdges),'rows'),'rows'),2);
            connM = sparse(edgesNotOuter(:,1),edgesNotOuter(:,2),ones(size(edgesNotOuter,1),1),max(edgesNotOuter(:)),max(edgesNotOuter(:)));
            connM = connM+connM';
            findConnector=@(x)C{x}(find(ismember(getIdsFromCC(C{x}),nonull(ff.segIds{ff_idx})),1));
            [~, path] = graphshortestpath(connM, findConnector(whichones(1)), findConnector(whichones(2)));
            connM(path,:)=0;
            connM(:,path)=0;
            theotherones = setdiff(1:4,whichones);
            [~, path2] = graphshortestpath(connM, C{theotherones(1)}(1), C{theotherones(2)}(1));
            rescued1 = [];
            rescued2 = [];
            filterNan0=@(x)x~=0&~isnan(x);
            if isempty(path2)
                rescued1 = find(filterNan0(graphshortestpath(connM, C{theotherones(1)}(1))));
                rescued2 = find(filterNan0(graphshortestpath(connM, C{theotherones(2)}(1))));
            end
            edgesToRescue = find(ismember(edgesNotOuter,sort([path(1:end-1)',path(2:end)';path2(1:end-1)',path2(2:end)'],2),'rows'));
                find(~ismember(sort([path(1:end-1)',path(2:end)';path2(1:end-1)',path2(2:end)'],2),edgesNotOuter,'rows'))
            assert(length(edgesToRescue)==size([path(1:end-1)',path(2:end)';path2(1:end-1)',path2(2:end)'],1));
            edgesToRescue = [edgesToRescue; find(all(ismember(edgesNotOuter,[rescued1,rescued2]),1))'];
            stillInLimbo = setdiff(setdiff(edgesNotOuter,lookupnodes(thisEdges),'rows'),edgesNotOuter(edgesToRescue,:),'rows');
            [~, probabilities] = ismember(stillInLimbo,graph.edges,'rows');
            probabilities(find(probabilities))=graph.prob(probabilities(find(probabilities)));
            [~, probidxs] = sort(probabilities,'descend'); %is that how it works
            counter = 1
            while true
                if counter > length(probidxs)
                    break;
                end
                
                flatten = @(x)x(:);
                if ~any(ismember(flatten([lookupnodes(thisEdges);edgesNotOuter(edgesToRescue,:)]), stillInLimbo(probidxs(counter),:)))
                    counter = counter + 1;
                    continue;
                end
                edgesTemp = [lookupnodes(thisEdges); edgesToRescue; stillInLimbo(probidxs(counter),:)];
                if length(Graph.findConnectedComponents(edgesTemp)) > 1
                    edgesToRescue=[edgesToRescue; stillInLimbo(probidxs(counter),:)];
                    
                else
                    edgesToDelete(end+1) = find(ismember(sort(thisEdgesCol,2),sort(stillInLimbo(probidxs(counter),:)),'rows'));
                end
                probidxs(counter)=[];   
            end
        end
        thisEdgesCol(edgesToDelete) = [];
        comments = cell(size([superagglo_definition.output.nodes;nodesToAdd],1), 1);
        writeNml([outputFolder 'skelSplit' num2str(idx) '.nml'], writeSkeletonFromNodesAndEdges({[superagglo_definition.output.nodes;nodesToAdd]}, {[thisEdgesCol;thisEdgesNew]}, {comments}, {'axon'}, {[0 0 1 1]}));
    end
end    
            
            
        
%        save(['/tmpscratch/kboerg/chiasmarunAugust/splitted_' num2str(idx)],'nodes2','edges2', 'segIds2');

function dummy()
end