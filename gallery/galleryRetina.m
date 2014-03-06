function galleryRetina()

addpath('/zdata/manuel/code/auxiliary/');
addpath('/zdata/manuel/code/auxiliary/hocMaker/');
addpath('/zdata/manuel/code/auxiliary/cubes/');

% Look up files
%files = dir('/zdata/manuel/sync/fromLap/ekSkel/bpc567subset/*.nml');
%load('/zdata/manuel/sync/fromLap/ekSkel/bpc567subset/skelsAsMat.mat');
files = dir('/zdata/manuel/sync/fromLap/ekSkel/bpcForPaperFinal/*.nml');

% Construct sphere for soma
%[x,y,z] = meshgrid(-55:55,-55:55,-55:55);
%se = (x/52).^2 + (y/52).^2 + (z/52).^2 <= 1;
%se = smooth3(se, 'gaussian', 5, 2);
%somaTemplate = isosurface(se);
%somaTemplate.vertices = somaTemplate.vertices - repmat(55,size(somaTemplate.vertices,1),3);
%somaTemplate.vertices = somaTemplate.vertices * 100;


for id=1:length(files)%length(skels)%length(files) %% CAREFUL: may be changed to do only part after error in nml
	%skel_data = {skels(id)};
	%nodes = skel_data{1}.nodes(:,1:3);
	skel_data = readNml(['/zdata/manuel/sync/fromLap/ekSkel/bpcForPaperFinal/' files(id).name]);
	nodes = skel_data{1,1}.nodes(:,1:3);
	if size(nodes,1) > 100
		% for each node, find cube in which it lays so the cubes data can be used for
		% several nodes
		nodeData.nodes = nodes;
		nodeData.cubeCoords = zeros(1,size(nodes,2));
		for i = 1 : size(nodes,1)
		   nodeData.cubeCoords(i,:) = floor(( nodes(i,:) - 1) / 128 ); %from readKnossosRoi ...  overlap of 128!
		   nodeData.flag(i) = 0;
		end
	
		%group nodes that lie in the same cube
		groupCount = 0;
		for i = 1 : size(nodeData.nodes,1)
		    if(nodeData.flag(i)==0)
		        groupCount = groupCount+1;
		        groupedNodes{1,groupCount}.nodes = nodeData.nodes(i,:);
		        groupedNodes{1,groupCount}.cubeCoords = nodeData.cubeCoords(i,:);
		        if(i<size(nodeData.nodes,1))            
		            for j = i+1 : size(nodeData.nodes,1)
		                if(nodeData.cubeCoords(i,:) == nodeData.cubeCoords(j,:))
		                    groupedNodes{1,groupCount}.nodes = [groupedNodes{1,groupCount}.nodes ; nodeData.nodes(j,:)];
		                    nodeData.flag(j) = 1;
		                end
		            end
		        end
		    end
		end
    
		% for each cube get data and write 1 into result cube for the voxels 
		% belonging to one of the segments of the nodes lying in the cube 
		for i = 1 : size(groupedNodes,2)
    			%read cube
			disp(['Skeleton ' num2str(id,'%.2i') ': processing cube ' num2str(i,'%.4i') ' of ' num2str(size(groupedNodes,2),'%.4i')]);
			%raw_cube = readKnossosCube('/zdata/manuel/data/cortex/2012-09-28_ex145_07x2/mag1/', '2012-09-28_ex145_07x2_mag1', groupedNodes{i}.cubeCoords, 'uint8', '', 'raw', 128);
			if all(groupedNodes{i}.cubeCoords > [3 2 0]) & all(groupedNodes{i}.cubeCoords < [31 40 43])
    			    cube = readKnossosCube('/nfs/bmo/mberning/mag1/', '100527_k0563_seg', groupedNodes{i}.cubeCoords, 'uint16', '', 'raw', 256);
			    %get the color values of the nodes in the cube
			    segIds = zeros(1,size(groupedNodes{i}.nodes,1));
			    zeroOfCube = groupedNodes{i}.cubeCoords * 128 + 1;
			    zeroOfCube = zeroOfCube-64;
			    for j = 1 : size(groupedNodes{i}.nodes,1) 
			        rel_coords = groupedNodes{i}.nodes(j,:) - zeroOfCube(1,:); % point in cube = actual point - zero of cube
			        segIds(1,j) = cube(rel_coords(1),rel_coords(2),rel_coords(3)); 
			    end   
			    segIds(segIds == 0) = []; %delete the zeros no neuron has color black
			    counts = histc(segIds,unique(segIds));
			    segIds = unique(segIds);
			    segIds(counts < 1) = [];
			    if ~isempty(segIds)
			    for k = 1 : size(segIds,2)
			        cube(cube == segIds(k)) = NaN;
			    end
			    cube(~isnan(cube)) = 0;
			    cube(isnan(cube)) = 1;
			    cube = imclose(cube, ones([3,3,3]));
			    cube = padarray(cube, [1 1 1]);
			    cube = smooth3(cube, 'gaussian', 5, 2);
			    issfs{i} = isosurface(cube, .5);
			    if ~isempty(issfs{i}.vertices)
				    reducepatch(issfs{i}, .2);
				    issfs{i}.vertices(:,[1 2]) = issfs{i}.vertices(:,[2 1]); 			    
				    issfs{i}.vertices = issfs{i}.vertices + repmat(zeroOfCube - [1 1 1],size(issfs{i}.vertices,1),1);
				    issfs{i}.vertices = issfs{i}.vertices .* repmat([12 12 25],size(issfs{i}.vertices,1),1);
			    end
			    end
        		end
		end
		if exist('issfs', 'var')
			idx = zeros(length(issfs),1);
			for i=1:length(issfs)
				if isempty(issfs{i}) | isempty(issfs{i}.vertices)
					idx(i) = 1;
				end
			end
			issfs(find(idx)) = [];
   			% ADD soma dirty version start
   			%[~,idx] = min(skel_data{1}.nodesNumDataAll(:,8));
    			%somaPos = skel_data{1}.nodesNumDataAll(idx,3:5);
    			%thisSoma = somaTemplate;
    			%thisSoma.vertices = somaTemplate.vertices + repmat(somaPos.*[12 12 25], size(somaTemplate.vertices,1), 1);
    			%issfs{end+1} = thisSoma;		
			% Save
			skel_data{1}.nodes(:,1:3) = skel_data{1}.nodes(:,1:3) .* repmat([12 12 25], size(skel_data{1}.nodes,1), 1);
			exportSurfaceToAmira(issfs, ['/zdata/manuel/sync/wholeCell/retina/bpc/' num2str(id, '%.2i') '.issf']);
			convertKnossosNmlToHoc(skel_data, ['/zdata/manuel/sync/wholeCell/retina/bpc/' num2str(id, '%.2i')], 0, 1, 0, 0, [1 1 1]);
			clear groupedNodes issfs;
		end
	end
end

end
