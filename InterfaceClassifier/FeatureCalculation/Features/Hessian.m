function imfeat = Hessian( raw, sigma, filterSize )
%HESSIAN Calculate the eigenvalues of the hessian of a 3D volume (image stack).
% INPUT sigma: Array containing standard deviation used for gaussian filter
%              in each dimension.
%       filterSize: Size of the resulting filter in each dimension and
%                   direction. The resulting filter has a size of
%                   2*filterSize + 1 and a total boundary of 2*filterSize.
% OUTPUT imfeat: cell array of length 3 containing the eigenvalues of the
%                hessian in increasing order (i.e. imfeat{3} is
%                the largest eigenvalue).
%
% NOTE There is a bug in the eigenvalue calculation (real and sort should
%      be interchanged). Use Codat.Features instead.
%
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

sigma = num2cell(sigma);
coords = arrayfun(@(siz)(-siz:siz).',filterSize,'UniformOutput',false);
coords = cellfun(@(coords,dim)permute(coords,[2:dim,1,(dim+1):3]),coords,num2cell(1:3),'UniformOutput',false);
gauss = cellfun(@(coords,sigma)exp(-(coords./sigma).^2./2),coords,sigma,'UniformOutput',false);
gauss = cellfun(@(gauss)gauss./sum(gauss),gauss,'UniformOutput',false);
gaussD = cellfun(@(coords,gauss,sigma)-gauss.*coords./(sigma.^2),coords,gauss,sigma,'UniformOutput',false);
gaussD2 = cellfun(@(coords,gauss,sigma)gauss.*(coords.^2./(sigma.^2)-1)./(sigma.^2),coords,gauss,sigma,'UniformOutput',false);
H = cell(3,3);
for dim1 = 1:3
    for dim2 = 1:dim1
        if(dim1==dim2)
            H{dim1,dim2} = convn(raw,gaussD2{dim1},'same');
        else
            H{dim1,dim2} = convn(raw,gaussD{dim1},'same');
            H{dim1,dim2} = convn(H{dim1,dim2},gaussD{dim2},'same');
        end
        for dim = setdiff(1:3,[dim1,dim2])
            H{dim1,dim2} = convn(H{dim1,dim2},gauss{dim},'same');
        end
        if (dim1 ~= dim2)
            H{dim2,dim1} = H{dim1,dim2};
        end
    end
end
[nx,ny,nz] = size(raw);
H = cellfun(@(x)x(:),H,'UniformOutput',false);
H = cellfun(@(x)permute(x,[3,2,1]),H,'UniformOutput',false);
%H = cell2mat(H);
%eigenvalues = eig3(H);
%eigenvalues = sort(eigenvalues,1);
%eigenvalues = real(eigenvalues);

%reduce RAM usage by processing H in blocks of size fraction
eigenvalues = zeros(3,nx*ny*nz,'single');
fraction = 3e7;
nParts = floor(nx*ny*nz/fraction);
for part=1:nParts
    HP = cellfun(@(x)x(1,1,(part-1)*fraction + 1:part*fraction),H,'UniformOutput',false);
    HP = cell2mat(HP);
    eigenvalues(:,(part-1)*fraction + 1:part*fraction) = eig3(HP);
    eigenvalues(:,(part-1)*fraction + 1:part*fraction) = sort(eigenvalues(:,(part-1)*fraction + 1:part*fraction),1);
    eigenvalues(:,(part-1)*fraction + 1:part*fraction) = real(eigenvalues(:,(part-1)*fraction + 1:part*fraction));
end
HP = cellfun(@(x)x(1,1,nParts*fraction + 1:end),H,'UniformOutput',false);
HP = cell2mat(HP);
clear H
eigenvalues(:,nParts*fraction + 1:end) = eig3(HP);
eigenvalues(:,nParts*fraction + 1:end) = sort(eigenvalues(:,nParts*fraction + 1:end),1);
eigenvalues(:,nParts*fraction + 1:end) = real(eigenvalues(:,nParts*fraction + 1:end));
clear HP

%save eigenvalues, where imfeat{1} contains the smallest
%eigenvalues and imfeat{3} the largest
imfeat = cell(3,1);
imfeat{1} = reshape(eigenvalues(1,:),nx,ny,nz);
imfeat{2} = reshape(eigenvalues(2,:),nx,ny,nz);
imfeat{3} = reshape(eigenvalues(3,:),nx,ny,nz);


end

