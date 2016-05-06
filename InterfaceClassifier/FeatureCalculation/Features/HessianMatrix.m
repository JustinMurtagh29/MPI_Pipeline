function H = HessianMatrix( raw, sigma, filterSize )
%HESSIAN Calculate the eigenvalues of the hessian of 3D data.
% INPUT sigma: Array containing standard deviation used for gaussian filter
%              in each dimension.
%       filterSize: Size of the resulting filter in each dimension and
%                   direction. The resulting filter has a size of
%                   2*filterSize + 1 and a total boundary of 2*filterSize.
% OUTPUT imfeat: cell array of length 3 containing the eigenvalues of the
%                hessian in increasing order (i.e. imfeat{3} is
%                the largest eigenvalue).

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
H = cellfun(@(x)x(:),H,'UniformOutput',false);
H = cellfun(@(x)permute(x,[3,2,1]),H,'UniformOutput',false);
H = cell2mat(H);


end

