function corrCoM(p)
%Calculate the CoM for the correspondence 'edges' [edges with p=1]

files = dir([p.saveFolder, 'correspondences/*.mat']);

for m = 1:length(files)

	load([p.saveFolder, '/correspondences/', files(m).name])
	load([p.local(result.cubeCoords1(1), result.cubeCoords1(1),result.cubeCoords1(3)).saveFolder, 'seg.mat']); 
	seg1 = seg;
	load([p.local(result.cubeCoords2(1), result.cubeCoords2(2),result.cubeCoords2(3)).saveFolder, 'seg.mat']); 
	seg2 = seg;
	
	for n = 1:size(result.correspondences,1)
		
		obj1 = seg1 == result.correspondences(n,1);
		obj2 = seg2 == result.correspondences(n,2);
		prop1 = regionprops(obj1, 'Centroid');
		prop2 = regionprops(obj2, 'Centroid');		
		CoM1(n,:) = num2cell(prop1.Centroid);
		CoM2(n,:) = num2cell(prop2.Centroid);
		CoM1 = cell2mat(CoM1);
		CoM2 = cell2mat(CoM2);

		%correction for global 3d coordinates depending on cubeCoords
		CoM1(:,1) = CoM1(:,1) + (result.cubeCoords1(1)-1) * p.tileSize(1);
		CoM1(:,2) = CoM1(:,2) + (result.cubeCoords1(2)-1) * p.tileSize(2);
		CoM1(:,3) = CoM1(:,3) + (result.cubeCoords1(3)-1) * p.tileSize(3);
		CoM2(:,1) = CoM2(:,1) + (result.cubeCoords2(1)-1) * p.tileSize(1);
		CoM2(:,2) = CoM2(:,2) + (result.cubeCoords2(2)-1) * p.tileSize(2);
		CoM2(:,3) = CoM2(:,3) + (result.cubeCoords2(3)-1) * p.tileSize(3);
		
	end

end


