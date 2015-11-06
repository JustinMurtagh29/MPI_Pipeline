function filtered = gaussianbpf(I,d0,d1)

% Fast n-dimensional fourier transform
[nx, ny, nz] = size(I);
fftI = fftn(I,[2*nx-1 2*ny-1 2*nz-1]);
fftI = fftshift(fftI);

% Initialize filter
[x,y,z] = ndgrid(-nx+1:nx-1,-ny+1:ny-1,-nz+1:nz-1);
dist = sqrt(x.^2 + y.^2 + z.^2);
clear x y z;
highPass = exp(-dist.^2./(2*d0^2));
lowPass = exp(-dist.^2./(2*d1^2));
clear dist;
bandPass = lowPass.*(1-highPass);
clear lowPass highPass;

% Update cube
filtered = bandPass.*fftI;
clear bandPass;
filtered = ifftshift(filtered);
filtered = ifftn(filtered,[2*nx-1 2*ny-1 2*nz-1]);
filtered = real(filtered(1:nx,1:ny,1:nz));

end
