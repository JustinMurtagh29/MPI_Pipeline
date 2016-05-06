function A = cropActivation( A, targetSize )
%CROPACTIVATION Crop activation to target size.
% INPUT A: 4D numerical array where the first 3 dimensions are cropped.
%       targetSize: [Nx1] array of integer specifying the target size of
%           the first three dimensions of A. The crop is made in the center
%           of A. If targetSize(4) > size(A,4) then the output is padded with
%           zeros in the fourth dimension to get the required target size
%           ("post" padding). targetSize(4) < size(A,4) is currently not
%           supported.
% OUTPUT A: The cropped input array.
%
% see also cnet.padArray
%
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

b = Codat.CNN.cnn.mSize(A,1:3) - targetSize(1:3);
b = floor(b./2);
A = A(b(1) + 1:b(1) + targetSize(1), ...
      b(2) + 1:b(2) + targetSize(2), ...
      b(3) + 1:b(3) + targetSize(3),:);
expandFMs = targetSize(4) - size(A,4);
if expandFMs > 0
    A = padarray(A,[0,0,0,expandFMs],'post');
end

end
