function calculateFeaturesCube(p, cubeNo, fm, voxelMap, filename)
% Calculate SynEM features for SegEM classification output on a local cube

% load segmentation
seg = loadSegDataGlobal(p.seg, p.local(cubeNo).bboxSmall);
if ~all(seg(:)==0)
   %load svg
   load(p.local(cubeNo).edgeFile, 'edges');
   load(p.local(cubeNo).borderFile, 'borders');
   % start empty
   features = [];
   if ~isempty(edges) && ~isempty(borders) && any([borders(:).Area] > fm.areaT)
	   %calculate interfaces
	   interfaces = SynEM.Svg.calculateInterfaces(seg, edges, borders, ...
		   fm.areaT, p.raw.voxelSize, fm.subvolsSize);
	   
	   %load voxel map
	   bboxFM = bsxfun(@plus, p.local(cubeNo).bboxSmall,[-fm.border', fm.border']./2);
	   switch lower(voxelMap)
		case 'class'
			vMap = loadClassData(p.class, bboxFM);
		case 'raw'
			vMap = loadRawData(p.raw, bboxFM);
            
            % NOTE(amotta): KMB's L2/3 dataset is significantly brighter
            % than his L4 dataset (by ~43 uint8). To train the classifiers
            % for L2/3, we fake the same distribution in L4.
            warning('Hard-coded offset for L2/3 simulation in L4');
            vMap = vMap + 43;
		case 'svm'
			vMap = loadSvmData(p, bboxFM);
		otherwise
			error('Unknown voxel map');
	   end
	   

	   % NOTE(amotta): It's possible for the voxel map to be four-dimensional.
		% This is, for example, the case for the output of Benedikt's SVM net.
		% We apply the same features map to all channels individually.

		for curIdx = 1:size(vMap, 4)
			%calculate features
			curFeatures = fm.calculate(interfaces, vMap(:, :, :, curIdx));

			% Save features
			if strcmp(fm.mode, 'direction')
				curFeatures = curFeatures(1:end/2,:); %only save first direction
			end
			features = cat(2, features, curFeatures);
		end
   end
   if ~exist('filename', 'var') || isempty(filename)
       filename = ['Interface' voxelMap 'Features.mat'];
   end
   outputFile = fullfile(p.local(cubeNo).saveFolder, filename);
   Util.save(outputFile, features);
end
end

function data = loadSvmData(param, box)
    % NOTE(amotta): This was copy-pasted from +Mito/loadPredictions of my
    % repository. Let's try to stay synchronized...
    
    % NOTE(amotta): At the moment, this is specific to the dataset
    % 2012-09-28_ex145_07x2_ROI2017. Once Benedikt's SVM CNN is a standard
    % part of the pipeline, this should be moved to the auxiliary methods.
    
    wkRoot = fullfile(param.saveFolder, 'svm');
    wkRoot = strcat(wkRoot, filesep);
    
    wkPrefix = '2012-09-28_ex145_07x2_corrected_mag1';
    wkCubeSize = [128, 128, 128, 4];
    
    data = readKnossosRoi( ...
        wkRoot, wkPrefix, box, 'single', '', 'raw', wkCubeSize);
end
