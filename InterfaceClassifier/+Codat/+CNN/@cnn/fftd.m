function out = fftd( A, siz, dim, d )
%FFTD Perform fft with d-regular array.
% INPUT A: Input array.
%       siz: Size of target in dim. Size is adapted to be processed
%            effectively by fft so check on output size.
%       dim: Apply fft along dim.
%       d: D-regularity factor.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

%adapt size to d-regularity factor
% r = rem(siz,d);
% if r > 0
%     siz = siz + d - r;
% end
siz = siz + mod(-siz,d);
reSize = ones(1,ndims(A));
reSize(dim) = d;
tmp = fft(A,siz/d,dim);
out = repmat(tmp,reSize);

end

