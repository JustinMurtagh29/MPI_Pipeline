function imfeat = SortedEigenvaluesHessian( raw, sigma, filterSize )
%SORTEDEIGENVALUESHESSIAN Calculate the eigenvalues of the hessian of a 3D volume (image stack).
% INPUT sigma: Array containing standard deviation used for gaussian filter
%              in each dimension.
%       filterSize: Size of the resulting filter in each dimension and
%                   direction. The resulting filter has a size of
%                   2*filterSize + 1 and a total boundary of 2*filterSize.
% OUTPUT imfeat: cell array of length 3 containing the eigenvalues of the
%                hessian in increasing order (i.e. imfeat{3} is
%                the larget eigenvalue).

sigma = num2cell(sigma);
coords = arrayfun(@(siz)(-siz:siz).',filterSize,'UniformOutput',false);
coords = cellfun(@(coords,dim)permute(coords,[2:dim,1,(dim+1):3]),coords,num2cell(1:3),'UniformOutput',false);
gauss = cellfun(@(coords,sigma)exp(-(coords./sigma).^2./2),coords,sigma,'UniformOutput',false);
gauss = cellfun(@(gauss)gauss./sum(gauss),gauss,'UniformOutput',false);
gaussD = cellfun(@(coords,gauss,sigma)-gauss.*coords./(sigma.^2),coords,gauss,sigma,'UniformOutput',false);
gaussD2 = cellfun(@(coords,gauss,sigma)gauss.*(coords.^2./(sigma.^2)-1)./(sigma.^2),coords,gauss,sigma,'UniformOutput',false);

Ixx1 = convn(raw,gaussD2{1},'same');
Ixx2 = convn(Ixx1,gauss{2},'same');
Ixx = convn(Ixx2,gauss{3},'same');
clear Ixx1 Ixx2

Iyy1 = convn(raw,gaussD2{2},'same');
Iyy2 = convn(Iyy1,gauss{1},'same');
Iyy = convn(Iyy2,gauss{3},'same');
clear Iyy1 Iyy2

Izz1 = convn(raw,gaussD2{3},'same');
Izz2 = convn(Izz1,gauss{1},'same');
Izz = convn(Izz2,gauss{2},'same');
clear Izz1 Izz2

Ixy1 = convn(raw,gaussD{1},'same');
Ixy2 = convn(Ixy1,gaussD{2},'same');
Ixy = convn(Ixy2,gauss{3},'same');
clear Ixy1 Ixy2

Ixz1 = convn(raw,gaussD{1},'same');
Ixz2 = convn(Ixz1,gaussD{3},'same');
Ixz = convn(Ixz2,gauss{2},'same');
clear Ixz1 Ixz2

Iyz1 = convn(raw,gaussD{2},'same');
Iyz2 = convn(Iyz1,gaussD{3},'same');
Iyz = convn(Iyz2,gauss{1},'same');
clear Ixy1 Ixy2

[nx,ny,nz] = size(raw);
imfeat = cell(1,3);
newSize = [1 1 nx*ny*nz];
eigenvalues = eig3([reshape(Ixx, newSize) reshape(Ixy, newSize) reshape(Ixz, newSize); ...
        reshape(Ixy, newSize) reshape(Iyy, newSize) reshape(Iyz, newSize); ...
        reshape(Ixz, newSize) reshape(Iyz, newSize) reshape(Izz, newSize)]);
eigenvalues = sort(eigenvalues,1);
eigenvalues = real(eigenvalues);
imfeat{1} = reshape(eigenvalues(1,:),nx,ny,nz);
imfeat{2} = reshape(eigenvalues(2,:),nx,ny,nz);
imfeat{3} = reshape(eigenvalues(3,:),nx,ny,nz);

end

