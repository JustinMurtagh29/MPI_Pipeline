function kl_cube = readKnossosCube( kl_parfolder, kl_fileprefix, kl_cubeCoord, classT,cubesize)

% READKNOSSOSCUBE: Read raw data from EM into Matlab
%
%   The function has the following arguments:
%       KL_PARFOLDER: Give the root directory of the data you want to read as a
%           string, e.g. 'E:\e_k0563\k0563_mag1\'
%       KL_FILEPREFIX: Give the name of the specific file you want to read without
%           the coordinates or the ending as a string, e.g. '100527_k0563_mag1'
%       KL_CUBECOORD: Give an array of 3 numbers of the xyz-coordinates as
%           given in the file name, no need for the full four digits: 0020 ->
%           20. E.g. [21 30 150]
%       CLASST: Optional! Standard version is unsigned int with 8 bits. For the
%           precision of the values.
%
%   => readKnossosCube( ?E:\e_k0563\k_0563_mag1', ?100527_k0563_mag1?, [21 30 150], ?uint8? )
%

if ~exist('classT','var') || isempty(classT)
    classT = 'uint8';
end

if ~exist('cubesize','var') || isempty(cubesize)
    cubesize=[128 128 128];
end
if numel(cubesize)==1
    cubesize=repmat(cubesize,1,3);
end

% Building the full filename
kl_fullfile = fullfile( kl_parfolder, sprintf( 'x%04.0f', kl_cubeCoord(1) ),...
    sprintf( 'y%04.0f', kl_cubeCoord(2) ), sprintf( 'z%04.0f', kl_cubeCoord(3) ),...
    sprintf( '%s_x%04.0f_y%04.0f_z%04.0f.raw', kl_fileprefix, kl_cubeCoord ) );

% If this file exists, load it into Matlab, else fill the matrix with zeros
if( exist( kl_fullfile, 'file' ) )
    fid = fopen( kl_fullfile,'r' );
    kl_cube = reshape(fread( fid,prod(cubesize), [classT '=>' classT]),cubesize);
    fclose( fid );
else
    kl_cube = zeros( cubesize, classT );
end
end