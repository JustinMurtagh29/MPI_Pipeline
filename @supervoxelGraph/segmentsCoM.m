function segCoM(p,i,j,k)
% calculates Center of Mass for unique segments connected by edges
% the output has numel(unique(edges(:))) which gives list with numRow = GlobalId of segment

        load(p.local(i,j,k).edgeFile)
	load(p.local(i,j,k).segFile)

	[uni, ia, ic] = unique([edges(:,1); edges(:,2)]);

%Calculate CoM of both supervoxels connected by edges
for m = 1:size(uni,1)

	object = seg == uni(m,1);   

	prop1 = regionprops(object, {'Centroid'});
	CoM(m,:) = num2cell(prop1.Centroid);
end

save([p.local(i,j,k).saveFolder, 'CoM.mat'], 'CoM');

end
