function [ deltaB, gW, gb ] = convLayerBwd( cnet, activity, delta, lyr, prop_down )
%CONVLAYERBWD Backward pass through convolutional layer.
% INPUT activity: Activity of layer lyr - 1
%       delta: Deltas (error derivative with respect to input) of layer
%       lyr.
%       prop_down: Bool indicating whether to backpropagate error to
%           previous layer. If false, only gradient wrt to parameters at
%           current layer are calculated.
%           (Optional. Default is false).
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if ~exist('prop_down','var') || isempty(prop_down)
    prop_down = false;
end

deltaB = zeros(size(activity),'like',activity);

if prod(cnet.filterSize{lyr}) == 1 %mlp layer
    %bias gradient
    gb = squeeze(sum(sum(sum(delta,1),2),3));
    %weight gradient
    if numel(activity)*size(delta,4)*4 < cnet.memoryLimit*1073741824 %use < 2GB memory
        gW = sum(sum(sum(bsxfun(@times,flip(activity,4),permute(delta,[1 2 3 5 4])),1),2),3) + ...
                          cnet.l2WeightDecayLambda.*cnet.W{lyr};
        if lyr > 2 %do not packpropagate to first layer
            deltaB = sum(bsxfun(@times,permute(delta,[1 2 3 5 4]),flip(cnet.W{lyr},4)),5);
        end
    else
        gW = zeros(size(cnet.W{lyr}),'like',cnet.W{lyr});
        for fm = 1:cnet.featureMaps(lyr)
            gW(:,:,:,:,fm) = cnet.flipdims(sum(sum(sum(bsxfun(@times,activity,delta(:,:,:,fm)),1),2),3)) + ...
                             cnet.l2WeightDecayLambda.*cnet.W{lyr}(:,:,:,:,fm);
            if lyr > 2
                deltaB = deltaB + bsxfun(@times,delta(:,:,:,fm),cnet.flipdims(cnet.W{lyr}(:,:,:,:,fm)));
            end
        end
    end
else %convolutional layer
    gb = shiftdim(sum(sum(sum(delta,1),2),3));
    switch cnet.convAlg
        case {'fft1','fft2'}
            %transform to fourier space
            actSize = size(activity);
            actSize(end+1:4) = 1;
            ffSize = actSize(1:3) + mod(-actSize(1:3),cnet.d(lyr,1:3));
            delSize = size(delta);
            delSizeNew = [ffSize(1:3), size(delta,4)];
            fftDimIdx = cnet.filterSize{lyr} > 1;
            ffSize(~fftDimIdx) = actSize(~fftDimIdx);
            fftDims = find(fftDimIdx);
            activity = cnet.flipdims(activity);
            for dim = fftDims
                activity = fft(activity,ffSize(dim),dim);
                delta = fft(delta,delSizeNew(dim),dim);
            end
            for dim = setdiff(1:3,fftDims)
                activity = flip(activity,dim);
            end
            gW = bsxfun(@times,activity,permute(delta,[1 2 3 5 4]));

            %gW calculation requires strided convolution with stride d
            %which should be implemented by calculating ifft with stride d
            %which is not natively implemented by matlab and hence full
            
            %ifft is calculated and unnecessary indices discarded)
            for dim = setdiff(1:3,fftDims) %accumulate gradient over batches
                gW = sum(gW,dim);
            end
            d = [cnet.d(lyr - 1,:), 1];
            fSize = cnet.filterSize{lyr};
            
            %fastest fft along first dim, indexing in other dim
            gW = ifft(gW,[],1); %1st dimension
            gW = gW(delSize(1):d(1):delSize(1) + d(1)*(fSize(1) - 1),:,:,:,:);
            gW = ifft(permute(gW,[2, 1, 3, 4, 5]),[],1); %second dimension
            if fftDimIdx(3)
                gW = permute(gW,[3, 2, 1, 4, 5]);
                gW = gW(:,:,delSize(2):d(2):delSize(2) + d(2)*(fSize(2) - 1),:,:);
                gW = ifft(gW,[],1);
                gW = permute(gW,[2, 3, 1, 4, 5]);
                gW = gW(:,:,delSize(3):d(3):delSize(3) + d(3)*(fSize(3) - 1),:,:);
            else
                gW = permute(gW,[2, 1, 3, 4, 5]);
                gW = gW(:,delSize(2):d(2):delSize(2) + d(2)*(fSize(2) - 1),:,:,:);
            end
            gW = real(gW) + cnet.l2WeightDecayLambda.*cnet.W{lyr};
            
            if prop_down
                W = flip(flip(flip(flip(cnet.W{lyr},1),2),3),4);
                delta = permute(delta,[1, 2, 3, 5, 4]);
                for dim = fftDims
                    W = cnet.fftd(W,ffSize(dim),dim,d(dim));
                end
                deltaB = sum(bsxfun(@times,delta,W),5);
                for dim = fftDims
                    deltaB = ifft(deltaB,ffSize(dim),dim);
                end
                deltaB = real(deltaB(1:actSize(1),1:actSize(2),1:actSize(3),:));
            end
            
        case 'convn'
            gW = zeros(size(cnet.W{lyr}),'like',cnet.W{lyr});
            for fm = 1:cnet.featureMaps(lyr)
                gW(:,:,:,:,fm) = convn(cnet.flipdims(activity),delta(:,:,:,fm),'valid') + ...
                                 cnet.l2WeightDecayLambda.*cnet.W{lyr}(:,:,:,:,fm);
                if prop_down
                    deltaB = deltaB + convn(delta(:,:,:,fm),cnet.flipdims(cnet.W{lyr}(:,:,:,:,fm)));
                end
            end
            %sparsify gradient (should be replaced by strided convolution
            %if speedup is necessary).
            gW = bsxfun(@times,gW,repmat(cnet.Wmask{lyr},[1, 1, 1, cnet.featureMaps(lyr - 1), cnet.featureMaps(lyr)]));
    end
end


end

