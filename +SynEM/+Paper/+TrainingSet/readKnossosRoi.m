function kl_roi = readKnossosRoi( kl_parfolder, kl_fileprefix, kl_bbox, classT,cubesize,overlap )
% READKNOSSOSROI: Read multiple raw data from EM into Matlab selecting a
%       region of interest
%
%   The function has the following arguments:
%       KL_PARFOLDER: Give the root directory of the data you want to read as a
%           string, e.g. 'E:\e_k0563\k0563_mag1\'
%       KL_FILEPREFIX: Give the name of the specific files you want to read without
%           the coordinates or the ending as a string, e.g. '100527_k0563_mag1'
%       KL_BBOX: Give an 3x2 array of the pixel coordinates of your region of
%           interest, no need for the full four digits: 0020 -> 20.
%           E.g. [129 384; 129 384; 129 384] CAREFUL: The amount of data will
%           easily explode.
%       CLASST: Optional! Standard version is unsigned int with 8 bits. For the
%           precision of the values.
%       CUBESIZE: Optional
%
%   => readKnossosRoi( ?E:\e_k0563\k0563_mag1', ?100527_k0563_mag1?,
%           [129 384; 129 384; 129 384], ?uint8? )
%
if ~exist('classT','var') || isempty(classT)
    classT = 'uint8';
end

if ~exist('cubesize','var') || isempty(cubesize)
    cubesize=128;
end
if numel(cubesize)==1
    cubesize=repmat(cubesize,1,3);
end
if ~exist('overlap','var') || isempty(overlap)
    overlap=0;
end
if numel(overlap)==1
    overlap=repmat(overlap,1,3);
end
ddims=cubesize(4:end);
% Calculate the size of the roi in pixels and xyz-coordinates
kl_bbox_size = kl_bbox(:,2)' - kl_bbox(:,1)' + [1 1 1];
kl_bbox_cubeind = [floor(( kl_bbox(:,1) - 1) / 128 ) ceil( kl_bbox(:,2) / 128 ) - 1];

kl_roi = zeros( [kl_bbox_size ddims], classT );
ncubes =diff(kl_bbox_cubeind,1,2)+1;
% Read every cube touched with readKnossosCube and write it in the right
% place of the kl_roi matrix
for i=1:prod(ncubes)
    kl_c=cell(3,1);
    [kl_c{:}]=ind2sub(ncubes.',i);
    kl_c=kl_bbox_cubeind(:,1)+cell2mat(kl_c)-1;
    
    kl_thiscube_coords = bsxfun(@times,[kl_c,kl_c+ 1] , cubesize.')+bsxfun(@times,overlap.',[-1,+1]);
    kl_thiscube_coords(:,1) = kl_thiscube_coords(:,1) + 1;
    
    kl_validbbox = [max( kl_thiscube_coords(:,1), kl_bbox(:,1) ),...
        min( kl_thiscube_coords(:,2), kl_bbox(:,2) )];
    
    kl_validbbox_cube = kl_validbbox - repmat( kl_thiscube_coords(:,1), [1 2] ) + 1;
    kl_validbbox_cubeC=cellfun(@(c)c(1):c(2),num2cell(kl_validbbox_cube,2),'UniformOutput',false);
    
    kl_validbbox_roi = kl_validbbox - repmat( kl_bbox(:,1), [1 2] ) + 1;
    kl_validbbox_roiC=cellfun(@(c)c(1):c(2),num2cell(kl_validbbox_roi,2),'UniformOutput',false);
    
    kl_cube = SynEM.Paper.TrainingSet.readKnossosCube( kl_parfolder, kl_fileprefix, kl_c, classT,cubesize+2*overlap );
    if(any(overlap))
        kl_roi( kl_validbbox_roiC{:},: ) = kl_roi( kl_validbbox_roiC{:},: )|(~kl_cube( kl_validbbox_cubeC{:},: ));
    else
        kl_roi( kl_validbbox_roiC{:},: ) = kl_cube( kl_validbbox_cubeC{:},: );
    end
end
 if(any(overlap))
    kl_roi=~kl_roi;
 end
end
