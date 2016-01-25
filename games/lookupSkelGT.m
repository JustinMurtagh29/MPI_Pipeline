function segIds = lookupSkelGT(p, skelGTflat)

    nodes = cat(1,skelGTflat.nodes{:});
    segIds = zeros(size(nodes,1),1, 'uint32');
    display('Reading segmentation IDs from KNOSSOS hierachy.');
    warning off; % To supress readKnossosRoi warning about empty cubes, should be uncommented/fixed otherwise
    tic;
    for x=1:size(p.local,1)
        for y=1:size(p.local,2)
            for z=1:size(p.local,3)
                % Determine which nodes are within small bounding box
                thisBbox = p.local(x,y,z).bboxSmall;
                nodeIdxInBbox = all(bsxfun(@ge, nodes, thisBbox(:,1)') & bsxfun(@le, nodes, thisBbox(:,2)'),2);
                % Load gloablized segmentation (and cut out small, non-overlapping part)
                load([p.local(x,y,z).saveFolder 'segGlobal.mat']);
                seg = seg(257:end-256,257:end-256,129:end-128);
                % Calculate position of nodes
                localPos = bsxfun(@minus, nodes(nodeIdxInBbox,:), thisBbox(:,1)' - [1 1 1]);
                if sum(segIds(nodeIdxInBbox)) ~= 0
                    error('WTF');
                end
                linearIdx = sub2ind(size(seg), localPos(:,1), localPos(:,2), localPos(:,3));
                segIds(nodeIdxInBbox) = seg(linearIdx);
            end
        end
        Util.progressBar(x,size(p.local,1));
    end
    warning on;
    segIds = mat2cell(segIds, cellfun(@(x)size(x,1), skelGTflat.nodes));

end
