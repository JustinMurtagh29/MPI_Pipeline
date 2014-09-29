classdef supervoxelGraph
	properties
		supervoxel
		supervoxelContacts
        adjMatrix
    end
    methods
		function sg = supervoxelGraph(p)
			sg.supervoxel = struct('cubeCoords', {}, 'segmentID', {}, 'CoM', {});
			sg.supervoxelContacts = struct('propability', {}, 'supervoxel1', {}, 'supervoxel2', {}, 'borderCoM', {});
            for i=1:size(p.local,1)
				for j=1:size(p.local,2)
					for k=1:size(p.local,3)
						load(p.local(i,j,k).edgeFile)
						load([p.local(i,j,k).saveFolder, 'CoM.mat'])
						load(p.local(i,j,k).borderFile)
                        load(p.local(i,j,k).probFile)

						idx = size(sg.supervoxel,2);
						[uni, ia, ic] = unique([edges(:,1); edges(:,2)]);						

						coords = cell(1,size(uni,1));
                        coords(1:end) = {[i,j,k]};
                        [sg.supervoxel(end+1:(end+size(coords,2))).cubeCoords] = coords{:};
                                                
                        segmentID = num2cell(uni(:,1))';
                        [sg.supervoxel((end-size(coords,2)+1) : end).segmentID] = segmentID{:};

						CoMcell = mat2cell(CoM, ones(size(CoM,1),1));   
						[sg.supervoxel((end-size(coords,2)+1) : end).CoM] = CoMcell{:};

						%create supervoxelContacts struct
						edgesGlobal = [ic(1:size(ic,1)/2), ic(size(ic,1)/2+1:end)] + idx;
						sv1 = num2cell(edgesGlobal(:,1));
						sv2 = num2cell(edgesGlobal(:,2));
						[sg.supervoxelContacts(end+1 : (end+size(sv1,1))).supervoxel1] = sv1{:};
						[sg.supervoxelContacts((end-size(sv1,1)+1) : end).supervoxel2] = sv2{:};

						bCoM = {borders(:).Centroid}';
                        [sg.supervoxelContacts((end-size(sv1,1)+1) : end).borderCoM] = bCoM{:};
						
						propC = num2cell(prob);
						[sg.supervoxelContacts((end-size(sv1,1)+1):end).propability] = propC{:};
			
					end
				end
			end

            %add edges from correspondences
%            edgesCorr = addCorrespondencesContacts(p, sg.supervoxel)
%            sv1 = num2cell(edgesCorr(:,1));
%            sv2 = num2cell(edgesCorr(:,2));
%            [sg.supervoxelContacts(end+1 : (end+size(sv1,1))).supervoxel1] = sv1{:};
%            [sg.supervoxelContacts((end-size(sv1,1)+1) : end).supervoxel2] = sv2{:};
%            [sg.supervoxelContacts((end-size(sv1,1)+1):end).propability] = 1;

            %create a NxN sparse AdjMatrix (best after addCorrespondenceEdges to have all supervoxel)
            sg.adjMatrix = sparse(size(sg.supervoxel,2), size(sg.supervoxel,2));
            idx = sub2ind(size(sg.adjMatrix), [sg.supervoxelContacts(:).supervoxel1], [sg.supervoxelContacts(:).supervoxel2]);
            idx2 = sub2ind(size(sg.adjMatrix), [sg.supervoxelContacts(:).supervoxel2], [sg.supervoxelContacts(:).supervoxel1]);
            sg.adjMatrix(idx) = [sg.supervoxelContacts.propability];
            sg.adjMatrix(idx2) = [sg.supervoxelContacts.propability];

            save([p.saveFolder, '/supervoxelGraph/supervoxel.mat'], 'sg')
        end

        function writeGlobalSegmentId(p)
        %write Segmentation with globalIds in knossos file structure and create array numEltotal for globalCorrId()
            globalSegID(p)      
        end

        function  globalCorrId(p)
        %create global mapping of corresponding segments w/ connComp
            
        end

%        function seededReconstruction(sg.adjMatrix, startSupervoxel)
%            output = seededReconstruction(sg.adjMatrix, startSupervoxel)
%        end
	end
end
