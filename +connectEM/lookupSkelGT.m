function [segIds, neighbours] = lookupSkelGT(p, skelGTflat)

    nodes = cat(1,skelGTflat.nodes{:});
    segIds = zeros(size(nodes, 1), 1, 'uint32');
    neighbours = zeros(size(nodes, 1), 26, 'uint32');
    tic;
    for x=1:size(p.local,1)
        for y=1:size(p.local,2)
            for z=1:size(p.local,3)
                % Determine which nodes are within small bounding box
                thisBbox = p.local(x,y,z).bboxSmall;
                nodeIdxInBbox = all(bsxfun(@ge, nodes, thisBbox(:,1)') & bsxfun(@le, nodes, thisBbox(:,2)'),2);
                if any(nodeIdxInBbox)
                    % Load gloablized segmentation
                    load([p.local(x,y,z).saveFolder 'segGlobal.mat']);
                    seg = padarray(seg, [1 1 1]);
                    % Calculate position of nodes
                    localPos = bsxfun(@minus, nodes(nodeIdxInBbox,:), thisBbox(:,1)' - [2 2 2]);
                    linearIdx = sub2ind(size(seg), localPos(:,1), localPos(:,2), localPos(:,3));
                    segIds(nodeIdxInBbox) = seg(linearIdx);
                    % Size of padded segmentation
                    [M,N,P] = size(seg);
                    % Construct 26-connectivity linear indices offset
                    neighIdxOffset = [(-M*N+[-M-1 -M -M+1 -1 0 1 M-1 M M+1]) [-M-1 -M -M+1 -1 1 M-1 M M+1] (M*N+[-M-1 -M -M+1 -1 0 1 M-1 M M+1])];
                    neighbours(nodeIdxInBbox,:) = seg(bsxfun(@plus, linearIdx, neighIdxOffset));
                end
            end
        end
        Util.progressBar(x, size(p.local,1));
    end
    nodesPerSkel = cellfun(@(x)size(x,1), skelGTflat.nodes);
    segIds = mat2cell(segIds, nodesPerSkel, 1);
    neighbours = mat2cell(neighbours, nodesPerSkel, 26);

end

