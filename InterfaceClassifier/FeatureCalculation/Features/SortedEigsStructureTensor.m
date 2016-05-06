function [ imfeat ] = SortedEigsStructureTensor( raw, sigmaW, filterSizeW, sigmaD, filterSizeD )
%SORTEDEIGSSTRUCTURETENSOR Calculate the eigenvalues of the Structure Tensor of 3D
% data.
% The total boundary of this filter is 2*filterSizeW + 2*filterSizeD
% INPUT sigmaW: Array containing standard deviation used for the window
%               size in each dimension.
%       filterSizeW: Size window in each dimension and direction.
%       sigmaD: Standard deviation used for the gradient in each dimension.
%       filterSizeD: Size of the gradient filter in each dimension and
%                    direction (see also GaussGradient).
% OUTPUT imfeat: cell array of length 3 containing the eigenvalues of the
%                structure tensor in increasing order (i.e. imfeat{3} is
%                the larget eigenvalue).

sigmaW = num2cell(sigmaW);
coords = arrayfun(@(siz)(-siz:siz).',filterSizeW,'UniformOutput',false);
coords = cellfun(@(coords,dim)permute(coords,[2:dim,1,(dim+1):3]),coords,num2cell(1:3),'UniformOutput',false);
gauss = cellfun(@(coords,sigma)exp(-(coords./sigma).^2./2),coords,sigmaW,'UniformOutput',false);
gauss = cellfun(@(gauss)gauss./sum(gauss),gauss,'UniformOutput',false);
grad = GaussGradient(raw,sigmaD,filterSizeD);

Ix = grad{1};
Iy = grad{2};
Iz = grad{3};
clear grad

Ixx1 = convn(Ix.^2,gauss{1},'same');
Ixx2 = convn(Ixx1,gauss{2},'same');
Ix_2 = convn(Ixx2,gauss{3},'same');
clear Ixx1 Ixx2

Iyy1 = convn(Iy.^2,gauss{1},'same');
Iyy2 = convn(Iyy1,gauss{2},'same');
Iy_2 = convn(Iyy2,gauss{3},'same');
clear Iyy1 Iyy2

Izz1 = convn(Iz.^2,gauss{1},'same');
Izz2 = convn(Izz1,gauss{2},'same');
Iz_2 = convn(Izz2,gauss{3},'same');
clear Izz1 Izz2

Ixy1 = convn(Ix.*Iy,gauss{1},'same');
Ixy2 = convn(Ixy1,gauss{2},'same');
Ixy = convn(Ixy2,gauss{3},'same');
clear Ixy1 Ixy2

Ixz1 = convn(Ix.*Iz,gauss{1},'same');
Ixz2 = convn(Ixz1,gauss{2},'same');
Ixz = convn(Ixz2,gauss{3},'same');
clear Ixz1 Ixz2

Iyz1 = convn(Iy.*Iz,gauss{1},'same');
Iyz2 = convn(Iyz1,gauss{2},'same');
Iyz = convn(Iyz2,gauss{3},'same');
clear Iyz1 Iyz2

[nx,ny,nz] = size(raw);
imfeat = cell(1,3);
newSize = [1 1 nx*ny*nz];
eigenvalues = eig3([reshape(Ix_2, newSize) reshape(Ixy, newSize) reshape(Ixz, newSize); ...
        reshape(Ixy, newSize) reshape(Iy_2, newSize) reshape(Iyz, newSize); ...
        reshape(Ixz, newSize) reshape(Iyz, newSize) reshape(Iz_2, newSize)]);
eigenvalues = sort(eigenvalues,1);
eigenvalues = real(eigenvalues);
imfeat{1} = reshape(eigenvalues(1,:),nx,ny,nz);
imfeat{2} = reshape(eigenvalues(2,:),nx,ny,nz);
imfeat{3} = reshape(eigenvalues(3,:),nx,ny,nz);

end

