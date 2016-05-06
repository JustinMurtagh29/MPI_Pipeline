function [ Y ] = maxPoolingBwd( cnet, X, ind, lyr )
%MAXPOOLINGBWD Backpropagation through 2x2x2 max pooling layer.
% This function readily allows for backpropagation through arbitrary
% pooling windows. This is not included yet because the max pooling forward
% pass is only through [2, 2, 2] pooling sizes.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

%calculate the global indices of the pooling window indices ind
d = cnet.d(lyr - 1,:).*cnet.stride{lyr};
if cnet.border(3) == 0
    d(3) = 0;
end
sX = size(X);
sX(end+1:4) = 1;
Y = zeros([sX(1:3) + d,sX(4)],'like',X);
sY = size(Y);
gl_ind = globalMpInd(ind, d, [2, 2, 2], sY);

%distribute the values in X to the acoording positions in Y
for k = 1:d(1) + 1
    for l = 1:d(2) + 1
        for m = 1:d(3) + 1
            Y(gl_ind(k:d(1) + 1:end,l:d(2) + 1:end,m:d(3)+1:end,:)) = ... 
                Y(gl_ind(k:d(1) + 1:end,l:d(2) + 1:end,m:d(3)+1:end,:)) +  ...
                X(k:d(1) + 1:end,l:d(2) + 1:end,m:d(3)+1:end,:);
        end
    end
end
end

function gl_ind = globalMpInd(mpInd, d, mpSize, sY)
% Calculate the global max-pooling indices

%indConv is used to include the sparsity of the mp window into the pooling
%indices
indConv = zeros(mpSize);
v = permute(0:mpSize(1) - 1,[2, 1]);
indConv = bsxfun(@plus,indConv,v.*d(1));
v = 0:mpSize(2) - 1;
indConv = bsxfun(@plus,indConv,v.*d(2)*sY(1));
v = permute(0:mpSize(3) - 1, [1 3 2]);
indConv = bsxfun(@plus,indConv,v.*d(3)*sY(1)*sY(2));
indConv = indConv(:);

%gl_indices use the linear first coordinate of each cube and convInd to
%convert the mpInd to linear indices in Y
poolingWindowInd = reshape(1:prod(sY),sY);
gl_ind = indConv(mpInd) + poolingWindowInd(1:end-d(1),1:end-d(2),1:end-d(3),:);

end

