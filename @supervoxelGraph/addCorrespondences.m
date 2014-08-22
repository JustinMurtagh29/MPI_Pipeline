function edgesCorr = addCorrespondences(p, supervoxel)
%Calculate the CoM for the correspondence 'edges' [edges with p=1]

files = dir([p.saveFolder, 'correspondences/*.mat']);

edgesCorr = [];

for m = 1:length(files)

	load([p.saveFolder, '/correspondences/', files(m).name])
	
	svCell = struct2cell(supervoxel);

	cubeCoords1 = [str2num(files(m).name(2)),str2num(files(m).name(4)), str2num(files(m).name(6))] 
	cubeCoords2 = [str2num(files(m).name(8)),str2num(files(m).name(10)), str2num(files(m).name(12))]; 
 	
	for n = 1:size(result.correspondences,1)
		globalID1(n) = find(cell2mat(arrayfun(@(x)(x.segmentID == result.correspondences(n,1) & all(x.cubeCoords == cubeCoords1) ),supervoxel,'uniformoutput',false)), 1, 'first'); 
		globalID2(n) = find(cell2mat(arrayfun(@(x)(x.segmentID == result.correspondences(n,2) & all(x.cubeCoords == cubeCoords2) ),supervoxel,'uniformoutput',false)), 1, 'first');
	end

	edgesCorr = [edgesCorr; globalID1', globalID2'];
end


