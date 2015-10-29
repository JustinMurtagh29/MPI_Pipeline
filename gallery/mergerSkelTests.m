function mergerSkelTests(p, skelPath, skelFile, outputPath)
    % Find probabilities of all edges merged together with this skeleton

    % Display status
    display(['Processing skeleton: ' skelFile]);
    % evalc to supress output of parseNml
    [~,skel_data] = evalc('parseNml([skelPath skelFile])');
    % Load supervoxel graph from disk 
    load([p.saveFolder 'graph.mat']);

    for sk = 1:length(skel_data)
        if size(skel_data{sk}.nodes,1) > 10
            nodeData.nodes = skel_data{sk}.nodes(:,1:3);
            nodeData.edges = skel_data{sk}.edges;
            pCube = p.local(:);
            bboxSmall = {pCube(:).bboxSmall};
            segIds = zeros(1,size(nodeData.nodes,1));
            % calculate local isosurfaces in global coordinates 
            for i=1:size(nodeData.nodes,1) 
                % Find ids of nodes
                segIds(i) = readKnossosRoi(p.seg.root, p.seg.prefix, [nodeData.nodes(i,:)+1' nodeData.nodes(i,:)'+1], 'uint32', '', 'raw'); 
            end
            for i=1:size(nodeData.edges,1)
                % segId1,2 for each node of current edge
                results(i).segId1 = segIds(nodeData.edges(i,1));
                results(i).segId2 = segIds(nodeData.edges(i,2));
                results(i).prob = graph.prob(any(graph.edges == results(i).segId1,2) & any(graph.edges == results(i).segId2,2));
            end
        end
    end
    save([outputPath strrep(skelFile, '.nml', '.mat')], 'results');
end

