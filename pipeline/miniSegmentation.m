function job = miniSegmentation(p)

for i=1:size(p.local,1)
	for j=1:size(p.local,2)
		for k=1:size(p.local,3)
			if ~exist(p.local(i,j,k).saveFolder, 'dir')
				mkdir(p.local(i,j,k).saveFolder);
			end
			idx = sub2ind(size(p.local), i, j, k);
			if isfield(p.local(i,j,k), 'class')
                % This is for pT train, probably needs update at some point?
				inputCell{idx} = {p.local(i,j,k).class.root, p.local(i,j,k).class.prefix, p.local(i,j,k).bboxBig, p.seg.func, p.local(i,j,k).segFile};
			else
				inputCell{idx} = {p.class.root p.class.prefix, p.local(i,j,k).bboxBig, p.seg.func, p.local(i,j,k).tempSegFile};
			end
		end
	end
end

functionH = @segmentForPipeline;
job = startCPU(functionH, inputCell, 'segmentation');

end

