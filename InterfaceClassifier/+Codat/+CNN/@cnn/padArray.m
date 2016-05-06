function A = padArray( A, targetSize )
%PADARRAY Pad array to target size.
% INPUT A: 4d array numerical array. The first three dimensions are padded.
%       targetSize: [Nx1] array of integer specifying the target size of
%           the first three dimensions of A. The original input array A is
%           at the center of the output. If targetSize(4) < size(A,4) then only
%           the first targetSize(4) indices of A along the fourth dimension
%           are used.
% OUTPUT A: The padded input array.
%
% see also cnet.cropActivation
%
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

sA = size(A);
b = targetSize(1:3) - sA(1:3);
b1 = floor(b./2);
if (size(A,4) - targetSize(4)) > 0
    A = A(:,:,:,1:targetSize(4));
end
A = padarray(A,[b1,0],'pre');
A = padarray(A,[b-b1,0],'post');

end
