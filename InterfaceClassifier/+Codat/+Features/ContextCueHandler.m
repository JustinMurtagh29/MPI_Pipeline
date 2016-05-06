classdef ContextCueHandler
    %CONTEXTCUEHANDLER Class for operating with context cues in 3d data.
    % Context cues are used in voxel-based processing. They provide
    % additional information around a voxel of interest by providing a set
    % of coordinates for each voxel of interest in a reference coordinate
    % system of that voxel. Context cues are thus placed consistently for
    % all voxel of interest and they can extend the feature space of a
    % a voxel by considering the feature responses also at the context cue
    % locations.
    %
    % PROPERTIES
    % numCues: The number of context cues.
    % maxCueNorm: The maximal distance of context cues to the reference
    %       voxel.
    % cues: The 3d locations of a context cue sample of size numCues and
    %       with maximal norm of maxCueNorm.
    %
    % Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>
    
    properties
        numCues
        maxCueNorm
        cues
    end
    
    methods
        function obj = ContextCueHandler(numCues, maxCueNorm)
            % Construct a context cue handler object.
            if nargin == 2
                obj.numCues = numCues;
                obj.maxCueNorm = maxCueNorm;
                obj = generateCueLocations(obj);
            else
                obj.numCues = [];
                obj.maxCueNorm = [];
                obj.cues = [];
            end
            
        end
        
        function obj = generateCueLocations(obj, numCues, maxCueNorm)
            %Generate context cues around the origin.
            % INPUT numCues: (Optional) The number of context cues
            %           locations to generate.
            %           (Default: obj.numCues if not empty)
            %       maxCueNorm: (Optional) The maximal allowed norm of each
            %           cue (i.e. the maximal distance to the point of
            %           interest).
            %           (Default: obj.maxCueNorm if not empty).
            
            if nargin == 1
                if isempty(obj.numCues) || isempty(obj.maxCueNorm)
                    error('NumCues or maxCueNorm is empty.');
                else
                    numCues = obj.numCues;
                    maxCueNorm = obj.maxCueNorm;
                end
            end
            [x,y,z] = meshgrid(-maxCueNorm:maxCueNorm, ...
                               -maxCueNorm:maxCueNorm, ...
                               -maxCueNorm:maxCueNorm);
            accessibleLocations = x.^2 + y.^2 + z.^2 <= maxCueNorm.^2;
            if sum(accessibleLocations(:)) < numCues
                error(['Maximal cue norm is too small to provide enough',...
                  'possible locations for the specified number of cues.']);
            end
            rndIdx = randperm(sum(accessibleLocations(:)),numCues);
            [x, y, z] = ind2sub(size(accessibleLocations), ...
                                find(accessibleLocations));
            cuesTmp = [x,y,z] - maxCueNorm - 1;
            cuesTmp = cuesTmp(rndIdx,:);
            obj.cues = cuesTmp;
        end
        
        function cues_li = calculateCueLocations(obj, raw, voxelList, type)
            % Calculate location of context cues for input voxels.
            % INPUT raw: 3d raw data array.
            %       voxelList: Array of linear indices of voxels of
            %           interest in raw.
            %       type: String specifying the calculation of the
            %           reference coordinate system. Options are
            %           'hessian','structureTensor'
            % OUTPUT cues_li: Matrix with length(voxelList) rows and
            %           obj.numCues column, where each row contains the
            %           linear index of the voxel of interest and the
            %           linear indices of all context cues for that voxel.
            
            if isrow(voxelList)
                voxelList = voxelList';
            end
            
            %get rotation matrix
            R = Codat.Features.ContextCueHandler.calculateCoordinateTransformation(raw, voxelList, type);
            
            %calculate cues locations
            [u,v,w] = ind2sub(size(raw),voxelList);
            coords = [u,v,w]';
            coords = mat2cell(coords,3,ones(length(voxelList),1))';
            refCues = [zeros(1,3); obj.cues];
            cues_li = cellfun(@(r,x)bsxfun(@plus,r*refCues',x),R,coords,'UniformOutput',false); %Becker formula (2)
            cues_li = cellfun(@(x)round(x),cues_li,'UniformOutput',false);
            cues_li = cellfun(@(x)sub2ind(size(raw),x(1,:),x(2,:),x(3,:)),cues_li,'UniformOutput',false);
            cues_li = cell2mat(cues_li);
        end
    end
    
    methods (Static)
        function out =getMatrixForVoxels(M, voxelList)
            out = zeros(length(voxelList),3,3,'like',M{1,1});
            for i = 1:3
                for j = i:3
                    out(:,i,j) = M{i,j}(voxelList);
                end
                for j = 1:i-1
                    out(:,i,j) = M{j,i}(voxelList);
                end
            end
        end
        
        function [R, D] = calculateCoordinateTransformation(raw, voxelList, type)
            %Calculate the rotation for the local coordinate transformation
            % around each voxel of interest.
            % INPUT raw: 3D raw data.
            %       voxelList: A vector of linear indices of the voxels of
            %           interest in raw.
            %       type: String specifying the calculation of the
            %           reference coordinate system. Options are
            %           'hessian','structureTensor'
            % OUTPUT R: Cell array where each cell contains the rotation
            %           matrix for the corresponding voxel of interest.
            %        D: Gauss gradient at the voxels of interest.
            
            if isrow(voxelList)
                voxelList = voxelList';
            end
            
            R = zeros(length(voxelList),3,3);
            raw = single(raw);
            switch type
                case 'hessian'
                    %sigma was chosen by inspecting some locations
                    sigma = 10./(2.*sqrt(2)).*[1 1 1]; %Becker below formula (11)
                    H = Codat.Features.hessianMatrix(raw,sigma);
                    H = Codat.Features.ContextCueHandler.getMatrixForVoxels(H, voxelList);
                    grad = Codat.Features.gaussGradient(raw,sigma);
                    D = [grad{1}(voxelList),grad{2}(voxelList),grad{3}(voxelList)];
                    [V,~] = eig3v3(H);
                    
                    %calculate rotation matrix (could be changed to bsxfun
                    %at some point)
                    for i = 1:length(voxelList)
                        R(i,:,3) = sign(V(i,:,3)*D(i,:)').*V(i,:,3); % Becker formula (9)
                        R(i,:,2) = sign(V(i,:,2)*D(i,:)').*V(i,:,2); % Becker formula (10)
                    end
                    R(:,:,1) = cross(R(:,:,3),R(:,:,2)); % Becker formula (11)
                    R = permute(R,[2 3 1]);
                    R = mat2cell(R,3,3,ones(length(voxelList),1));
                    R = squeeze(R);
                case 'structureTensor'
                    sigma = 10./(2.*sqrt(2)).*[1 1 1];
                    S = Codat.Features.structureTensor(raw ,sigma, sigma);
                    S = Codat.Features.ContextCueHandler.getMatrixForVoxels(S, voxelList);
                    grad = Codat.Features.gaussGradient(raw,sigma);
                    D = [grad{1}(voxelList),grad{2}(voxelList),grad{3}(voxelList)];
                    [V,~] = eig3v3(S);
                    
                    %calculate rotation matrix
                    for i = 1:length(voxelList)
                        R(i,:,3) = sign(V(i,:,3)*D(i,:)').*V(i,:,3);
                        R(i,:,2) = sign(V(i,:,2)*D(i,:)').*V(i,:,2);
                    end
                    R(:,:,1) = cross(R(:,:,3),R(:,:,2));
                    R = permute(R,[2 3 1]);
                    R = mat2cell(R,3,3,ones(length(voxelList),1));
                    R = squeeze(R);
                otherwise
                    error('Type %s not implemented.', type)
            end
        end
    end
    
end

