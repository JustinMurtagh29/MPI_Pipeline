function activityWithoutNL = convLayerFwd( cnet, activity, lyr )
%CONVLAYERFWD Forward pass through convolutional layer.
% INPUT activity: Activity of layer lyr-1
%       lyr: Current layer.
% OUTPUT activityWithoutNL: Output of convolutional layer.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

activityWithoutNL = zeros([cnet.mSize(activity,1:3) - (cnet.dFilterSize{lyr} - 1),cnet.featureMaps(lyr)],'like',activity);
%mlpConv layer (do direct multiplication instead of fft/convn)
if prod(cnet.filterSize{lyr}) == 1  && numel(activity)*size(cnet.W{lyr},5)*4 < cnet.memoryLimit*1073741824
    activityWithoutNL = bsxfun(@plus,squeeze(sum(bsxfun(@times,activity,flip(cnet.W{lyr},4)),4)),permute(cnet.b{lyr},[4 3 2 1]));
elseif prod(cnet.filterSize{lyr}) == 1
    for fm = 1:cnet.featureMaps(lyr)
        activityWithoutNL(:,:,:,fm) = sum(bsxfun(@times,activity,flip(cnet.W{lyr}(:,:,:,:,fm),4)),4) + cnet.b{lyr}(fm);
    end
else %convolutional layer
    switch cnet.convAlg
        case 'fft1' %fully vectorized and adapted for strided convolutions
            %valid convolution via fft over whole array (high memory
            %requirement)
            
            %transform to fourier space
            W = cnet.W{lyr}(:,:,:,end:-1:1,:); %compatibility with old version
            actSize = cnet.mSize(activity,1:4);
            ffSize = [actSize(1:3) + mod(-actSize(1:3),cnet.d(lyr,1:3)),actSize(4)];
            wSize = cnet.dFilterSize{lyr};
            d = [cnet.d(lyr - 1,:), 1];
            %only do fft if filter size > 1
            fftDims = cnet.filterSize{lyr} > 1;
            ffSize(~fftDims) = actSize(~fftDims);
            fftDims = find(fftDims);
            for k = fftDims
                activity = fft(activity,ffSize(k),k);
                W = cnet.fftd(W,ffSize(k),k, d(k));
            end
            
            %fully vectorized product in fourier space
%             C = bsxfun(@times,activity,W);
%             C = permute(sum(C,4),[1, 2, 3, 5, 4]);

            %semi-vectorized with usually much less memory requirement and
            %almost same speed
            C = zeros([ffSize(1:3), cnet.featureMaps(lyr)],'like',activity);
            for i = 1:cnet.featureMaps(lyr)
                C(:,:,:,i) = sum(bsxfun(@times,activity,W(:,:,:,:,i)),4); %bsxfun necessary if any dimension in w is 1
            end
            
            for k = fftDims
                C = ifft(C,ffSize(k),k);
            end
            activityWithoutNL = bsxfun(@plus,real(C(wSize(1):actSize(1),wSize(2):actSize(2),wSize(3):actSize(3),:)),permute(cnet.b{lyr},[4 3 2 1]));

        case 'fft2'
            %valid convolution via fft iterated over fms (medium memory
            %requirement)
            W = cnet.W{lyr}(:,:,:,end:-1:1,:); %compatibility with old version
            szA = cnet.mSize(activity,1:4);
            szW = cnet.dFilterSize{lyr};
            ffSize = [szA(1:3) + mod(-szA(1:3),cnet.d(lyr,1:3)),szA(4)];
            d = [cnet.d(lyr - 1,:), 1];
            
            %only do fft if filter size > 1
            fftDims = cnet.filterSize{lyr} > 1;
            ffSize(~fftDims) = szA(~fftDims);
            fftDims = find(fftDims);
            
            %fft activity
            for k = fftDims
                activity = fft(activity,ffSize(k),k);
            end
            
            %fft weights and ifft activity*weights for current feature map
            for fm = 1:cnet.featureMaps(lyr)
                w = W(:,:,:,:,fm);
                for k = fftDims
                    w = cnet.fftd(w,ffSize(k),k, d(k));
                end
                C = sum(bsxfun(@times,activity,w),4); %bsxfun necessary if any dimension in w is 1
                for k = fftDims
                    C = ifft(C,ffSize(k),k);
                end
                activityWithoutNL(:,:,:,fm) = real(C(szW(1):szA(1),szW(2):szA(2),szW(3):szA(3))) + cnet.b{lyr}(fm);
            end

        case 'convn'
            %standard convolution iterated over fms (lowest memory requirement)
            for fm = 1:cnet.featureMaps(lyr)
                activityWithoutNL(:,:,:,fm) = convn(activity,cnet.W{lyr}(:,:,:,:,fm),'valid') + cnet.b{lyr}(fm);
            end
        otherwise
            error('Unknown algorithm %s.',cnet.convAlg);
    end
end

end

