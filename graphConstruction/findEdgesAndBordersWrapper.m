function findEdgesAndBordersWrapper(segFile, edgeFile, borderFile, segmentFile, raw, bbox)
    % Computation of edges and borders optimized for 512x512x256

    % Load *SMALL* cube of global segmentation IDs
    load(segFile);

    % Use new function from SynEM for calculation of svg data
    if nargin > 4
	mask = readKnossosRoi(regexprep(raw.root,'/color/','/mask/'),[raw.prefix,'_mask'],bbox);
        if  all(mask(:))
		[edges, borders, segments] = SynEM.Svg.findEdgesAndBorders(seg);
	elseif ~any(mask(:))
		edges =  zeros(0,2,'uint32');
		borders = struct('PixelIdxList',{},'Area',{},'Centroid',{});
		segments = struct('PixelIdxList',{},'Id',{});
	else
		[edges, borders, segments] = SynEM.Svg.findEdgesAndBorders(seg,~mask);
	end
%        if all(excludeVoxels(:)) % skip calculation if completely in mirrorPaded region
%            return
%        end
%        [edges, borders, segments] = SynEM.Svg.findEdgesAndBorders(seg,excludeVoxels);
    else
        [edges, borders, segments] = SynEM.Svg.findEdgesAndBorders(seg);
    end

    % Save to files: currently segmentation overwrites old one (now leaves are merged)
    Util.save(edgeFile, edges);
    Util.save(borderFile, borders);
    Util.save(segmentFile, segments);
end
