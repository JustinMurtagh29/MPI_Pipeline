function calcFeatures(parameter, sub)

idx = sub2ind(sub,1);

if ~exist(parameter.local(sub(1),sub(2),sub(3)).saveFolder, 'dir') 
	mkdir(parameter.local(sub(1),sub(2),sub(3)).saveFolder);
end

if strcmp(parameter.feature.input, 'raw')
	imfeat = loadRawData(parameter.raw.root, parameter.raw.prefix, parameter.local(sub(1),sub(2),sub(3)).bboxSmall, 1)                
end

if strcmp(parameter.feature.input, 'aff')
	imfeat = loadClassData(parameter.class.root, parameter.class.prefix, parameter.local(sub(1),sub(2),sub(3)).bboxSmall);
end

border = load(parameter.local(sub(1),sub(2),sub(3)).borderFile);

for m=1:size(parameter.filter,2)
	imfeats = filter3d(parameter, imfeat, m);
	weights = featureDesign(imfeats, border);
	for l = 1:7
        	weights(:,(l-1)*7+l) = weights(:,l);                               %7 Features/Filter;   
	end
end

save([parameter.local(sub(1),sub(2),sub(3)).saveFolder  'weight-' num2str(idx) '-' num2str(i) '.mat'], 'weights');

end
