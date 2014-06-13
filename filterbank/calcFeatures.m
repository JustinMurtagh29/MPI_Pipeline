function calcFeatures(parameter, sub)

idx = sub2ind(sub,1);
weights = [];

if ~exist(parameter.local(sub(1),sub(2),sub(3)).saveFolder, 'dir') 
	mkdir(parameter.local(sub(1),sub(2),sub(3)).saveFolder);
end

for l=1:length(parameter.feature.input) 
	if strcmp(parameter.feature.input{l}, 'raw')
		imfeat = loadRawData(parameter.raw.root, parameter.raw.prefix, parameter.local(sub(1),sub(2),sub(3)).bboxSmall, 1);           
	elseif strcmp(parameter.feature.input{l}, 'aff')
		imfeat = loadClassData(parameter.class.root, parameter.class.prefix, parameter.local(sub(1),sub(2),sub(3)).bboxSmall);
	end
	load(parameter.local(sub(1),sub(2),sub(3)).borderFile);
	for m=1:size(parameter.filter,2)
		for n=1:length(parameter.filter{m}{2})
	 		imfeats = filter3d(parameter, imfeat, m, n);
 	 		if isa(imfeats, 'cell')
				for p=1:length(imfeats)
					weights_new = featureDesign(real(imfeats{p}), borders);
					weights = [weights weights_new];
				end
	       		else
				weights_new = featureDesign(imfeats, borders);
				weights = [weights weights_new];	
			end
		end
	end
end

save(parameter.local(sub(1),sub(2),sub(3)).weightFile, 'weights');

end

