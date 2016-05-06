function imfeat = StructureTensor( raw, sigmaW, filterSizeW, sigmaD, filterSizeD )
%STRUCTURETENSOR Calculate the eigenvalues of the Structure Tensor of 3D
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
%                the largest eigenvalue).
%
% NOTE There is a bug in the eigenvalue calculation (real and sort should
%      be interchanged). Use Codat.Features instead.
%
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

sigmaW = num2cell(sigmaW);
coords = arrayfun(@(siz)(-siz:siz).',filterSizeW,'UniformOutput',false);
coords = cellfun(@(coords,dim)permute(coords,[2:dim,1,(dim+1):3]),coords,num2cell(1:3),'UniformOutput',false);
gauss = cellfun(@(coords,sigma)exp(-(coords./sigma).^2./2),coords,sigmaW,'UniformOutput',false);
gauss = cellfun(@(gauss)gauss./sum(gauss),gauss,'UniformOutput',false);
grad = GaussGradient(raw,sigmaD,filterSizeD);
S = cell(3,3);
for dim1 = 1:3
    for dim2 = 1:dim1
        S{dim1,dim2} = grad{dim1}.*grad{dim2};
        for dim3=1:3
            S{dim1,dim2} = convn(S{dim1,dim2},gauss{dim3},'same');
        end
        if (dim1 ~= dim2)
            S{dim2,dim1} = S{dim1,dim2};
        end
    end
end
[nx,ny,nz] = size(raw);
S = cellfun(@(x)x(:),S,'UniformOutput',false);
S = cellfun(@(x)permute(x,[3,2,1]),S,'UniformOutput',false);
%S = cell2mat(S);
%eigenvalues = eig3(S);
%eigenvalues = sort(eigenvalues,1);
%eigenvalues = real(eigenvalues);

%reduce RAM usage by processing H in blocks of size fraction
eigenvalues = zeros(3,nx*ny*nz,'single');
fraction = 3e7;
nParts = floor(nx*ny*nz/fraction);
for part=1:nParts
    SP = cellfun(@(x)x(1,1,(part-1)*fraction + 1:part*fraction),S,'UniformOutput',false);
    SP = cell2mat(SP);
    eigenvalues(:,(part-1)*fraction + 1:part*fraction) = eig3(SP);
    eigenvalues(:,(part-1)*fraction + 1:part*fraction) = sort(eigenvalues(:,(part-1)*fraction + 1:part*fraction),1);
    eigenvalues(:,(part-1)*fraction + 1:part*fraction) = real(eigenvalues(:,(part-1)*fraction + 1:part*fraction));
end
SP = cellfun(@(x)x(1,1,nParts*fraction + 1:end),S,'UniformOutput',false);
SP = cell2mat(SP);
clear S
eigenvalues(:,nParts*fraction + 1:end) = eig3(SP);
eigenvalues(:,nParts*fraction + 1:end) = sort(eigenvalues(:,nParts*fraction + 1:end),1);
eigenvalues(:,nParts*fraction + 1:end) = real(eigenvalues(:,nParts*fraction + 1:end));
clear SP

%save eigenvalues, where imfeat{1} contains the smallest
%eigenvalues and imfeat{3} the largest
imfeat = cell(3,1);
imfeat{1} = reshape(eigenvalues(1,:),nx,ny,nz);
imfeat{2} = reshape(eigenvalues(2,:),nx,ny,nz);
imfeat{3} = reshape(eigenvalues(3,:),nx,ny,nz);

end

