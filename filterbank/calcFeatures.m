function calcFeatures(parameter, sub)

idx = sub2ind(sub,1);

if ~exist(parameter.local(sub(1),sub(2),sub(3)).saveFolder, 'dir') 
	mkdir(parameter.local(sub(1),sub(2),sub(3)).saveFolder);
end

if strcmp(parameter.feature.input, 'raw')
	imfeat = loadRawData(parameter.raw.root, parameter.raw.prefix, parameter.local(sub(1),sub(2),sub(3)).bboxBig, 1)                
end

if strcmp(parameter.feature.input, 'aff')
	imfeat = loadClassData(parameter.class.root, parameter.class.prefix, parameter.local(sub(1),sub(2),sub(3)).bboxBig);
end

border = load('~/GIT/manuelCode/border.mat');
border = border.border;
disp('border loaded');

for m=1:1%size(parameter.filter,2)
	imfeats = filter3d(parameter, imfeat, m);
    	disp('filter calculated');
	weights = featureDesign(imfeats, border);
	disp('weights calculated');
    for l = 1:7
        weights(:,(l-1)*7+l) = weights(:,l);                               %7 Features/Filter;   
    end
 disp(m);
end

save([parameter.local(sub(1),sub(2),sub(3)).saveFolder  'weight-' num2str(idx) '-' num2str(i) '.mat'], 'weights');
