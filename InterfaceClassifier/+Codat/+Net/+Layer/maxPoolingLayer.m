classdef maxPoolingLayer < Codat.Net.Layer.baseLayer
    %MAXPOOLINGLAYER Max pooling layer for dense pixel/voxel classification.
    % The max pooling layer calculates the maximum value from a pooling window
    % shifted over the input image.
    %
    % PROPERTIES
    % mpSize: Integer vector of length 3 defining the size of the pooling window
    %         in each dimension.
    % sparsityFactor: Integer vector of length 3 defining the sparsity of the
    %         pooling window. This is only important for dense pixel
    %         classification and depends on the strides of previous layers.
    %         The sparsity factor of this layer should be the product of
    %         stides of preceding layers.
    % stride: Integer vector of length 3 defining the stride in each dimension
    %         (i.e. the number of pixels of a single shift of the pooling
    %         window). For dense pixel classification this value should be set
    %         to [1, 1, 1] and the sparsityFactor of succeding layers should be
    %         multiplied by the desired stride.
    % mpInd: Linear indices of the field in the pooling window containing the
    %        maximum. If the pooling window is sparse, then only
    % Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>
    
    properties
        mpSize;
        sparsityFactor;
        mpInd;
    end
    
    methods
        function obj = maxPoolingLayer( mpSize, sparsityFactor, stride)
            obj.mpSize = mpSize;
            if ~exist('sparsityFactor','var') || isempty(sparsityFactor)
                obj.sparsityFactor = [1, 1, 1];
            else
                obj.sparsityFactor = sparsityFactor;
            end
            if ~exist('stride','var') || isempty(stride)
                obj.stride = [1, 1, 1];
            else
                obj.stride = stride;
            end
            obj.isInPlace = true;
        end
        
        function [ layer, Y ] = forward(layer, X)
            d = layer.sparsityFactor;
            
            %some frequently used max pooling windows are hardcoded which
            %is both faster and less memory consuming than the general
            %pooling version (which should be optimized at some point)
            if all(layer.mpSize == [2, 2, 2])
                [Y,ind] = max([subsref(X(1:end-d(1),1:end-d(2),1:end-d(3),:),struct('type','()','subs',{{':'}})), ...
                    subsref(X(1 + d(1):end,1:end-d(2),1:end-d(3),:),struct('type','()','subs',{{':'}})), ...
                    subsref(X(1:end-d(1),1 + d(2):end,1:end-d(3),:),struct('type','()','subs',{{':'}})), ...
                    subsref(X(1 + d(1):end,1 + d(2):end,1:end-d(3),:),struct('type','()','subs',{{':'}})), ...
                    subsref(X(1:end-d(1),1:end-d(2),1 + d(3):end,:),struct('type','()','subs',{{':'}})), ...
                    subsref(X(1 + d(1):end,1:end-d(2),1 + d(3):end,:),struct('type','()','subs',{{':'}})), ...
                    subsref(X(1:end-d(1),1 + d(2):end,1 + d(3):end,:),struct('type','()','subs',{{':'}})), ...
                    subsref(X(1 + d(1):end,1 + d(2):end,1 + d(3):end,:),struct('type','()','subs',{{':'}}))],[],2);
            elseif all(layer.mpSize == [2, 2, 1])
                [Y,ind] = max([subsref(X(1:end-d(1),1:end-d(2),:,:),struct('type','()','subs',{{':'}})), ...
                    subsref(X(1 + d(1):end,1:end-d(2),:,:),struct('type','()','subs',{{':'}})), ...
                    subsref(X(1:end-d(1),1 + d(2):end,:,:),struct('type','()','subs',{{':'}})), ...
                    subsref(X(1 + d(1):end,1 + d(2):end,:,:),struct('type','()','subs',{{':'}}))],[],2);
            elseif all(layer.mpSize == [3, 3, 3])
                d2 = 2.*d;
                [Y,ind] = max([subsref(X(1:end-d2(1),1:end-d2(2),1:end-d2(3),:),struct('type','()','subs',{{':'}})), ...
                    subsref(X(1 + d(1):end - d(1),1:end-d2(2),1:end-d2(3),:),struct('type','()','subs',{{':'}})), ...
                    subsref(X(1 + d2(1):end,1:end-d2(2),1:end-d2(3),:),struct('type','()','subs',{{':'}})), ...
                    subsref(X(1:end-d2(1),1 + d(2):end - d(2),1:end-d2(3),:),struct('type','()','subs',{{':'}})), ...
                    subsref(X(1 + d(1):end - d(1),1 + d(2):end - d(2),1:end-d2(3),:),struct('type','()','subs',{{':'}})), ...
                    subsref(X(1 + d2(1):end,1 + d(2):end - d(2),1:end-d2(3),:),struct('type','()','subs',{{':'}})), ...
                    subsref(X(1:end - d2(1),1 + d2(2):end,1:end-d2(3),:),struct('type','()','subs',{{':'}})), ...
                    subsref(X(1 + d(1):end - d(1),1 + d2(2):end,1:end-d2(3),:),struct('type','()','subs',{{':'}})), ...
                    subsref(X(1 + d2(1):end,1 + d2(2):end,1:end-d2(3),:),struct('type','()','subs',{{':'}})), ...
                    subsref(X(1:end-d2(1),1:end-d2(2),1 + d(3):end-d(3),:),struct('type','()','subs',{{':'}})), ...
                    subsref(X(1 + d(1):end - d(1),1:end-d2(2),1 + d(3):end-d(3),:),struct('type','()','subs',{{':'}})), ...
                    subsref(X(1 + d2(1):end,1:end-d2(2),1 + d(3):end-d(3),:),struct('type','()','subs',{{':'}})), ...
                    subsref(X(1:end-d2(1),1 + d(2):end - d(2),1 + d(3):end-d(3),:),struct('type','()','subs',{{':'}})), ...
                    subsref(X(1 + d(1):end - d(1),1 + d(2):end - d(2),1 + d(3):end-d(3),:),struct('type','()','subs',{{':'}})), ...
                    subsref(X(1 + d2(1):end,1 + d(2):end - d(2),1 + d(3):end-d(3),:),struct('type','()','subs',{{':'}})), ...
                    subsref(X(1:end - d2(1),1 + d2(2):end,1 + d(3):end-d(3),:),struct('type','()','subs',{{':'}})), ...
                    subsref(X(1 + d(1):end - d(1),1 + d2(2):end,1 + d(3):end-d(3),:),struct('type','()','subs',{{':'}})), ...
                    subsref(X(1 + d2(1):end,1 + d2(2):end,1 + d(3):end-d(3),:),struct('type','()','subs',{{':'}})), ...
                    subsref(X(1:end-d2(1),1:end-d2(2),1 + d2(3):end,:),struct('type','()','subs',{{':'}})), ...
                    subsref(X(1 + d(1):end - d(1),1:end-d2(2),1 + d2(3):end,:),struct('type','()','subs',{{':'}})), ...
                    subsref(X(1 + d2(1):end,1:end-d2(2),1 + d2(3):end,:),struct('type','()','subs',{{':'}})), ...
                    subsref(X(1:end-d2(1),1 + d(2):end - d(2),1 + d2(3):end,:),struct('type','()','subs',{{':'}})), ...
                    subsref(X(1 + d(1):end - d(1),1 + d(2):end - d(2),1 + d2(3):end,:),struct('type','()','subs',{{':'}})), ...
                    subsref(X(1 + d2(1):end,1 + d(2):end - d(2),1 + d2(3):end,:),struct('type','()','subs',{{':'}})), ...
                    subsref(X(1:end - d2(1),1 + d2(2):end,1 + d2(3):end,:),struct('type','()','subs',{{':'}})), ...
                    subsref(X(1 + d(1):end - d(1),1 + d2(2):end,1 + d2(3):end,:),struct('type','()','subs',{{':'}})), ...
                    subsref(X(1 + d2(1):end,1 + d2(2):end,1 + d2(3):end,:),struct('type','()','subs',{{':'}}))],[],2);
            elseif all(layer.mpSize == [3, 3, 2])
                d2 = 2.*d;
                [Y,ind] = max([subsref(X(1:end-d2(1),1:end-d2(2),1:end-d(3),:),struct('type','()','subs',{{':'}})), ...
                    subsref(X(1 + d(1):end - d(1),1:end-d2(2),1:end-d(3),:),struct('type','()','subs',{{':'}})), ...
                    subsref(X(1 + d2(1):end,1:end-d2(2),1:end-d(3),:),struct('type','()','subs',{{':'}})), ...
                    subsref(X(1:end-d2(1),1 + d(2):end - d(2),1:end-d(3),:),struct('type','()','subs',{{':'}})), ...
                    subsref(X(1 + d(1):end - d(1),1 + d(2):end - d(2),1:end-d(3),:),struct('type','()','subs',{{':'}})), ...
                    subsref(X(1 + d2(1):end,1 + d(2):end - d(2),1:end-d(3),:),struct('type','()','subs',{{':'}})), ...
                    subsref(X(1:end - d2(1),1 + d2(2):end,1:end-d(3),:),struct('type','()','subs',{{':'}})), ...
                    subsref(X(1 + d(1):end - d(1),1 + d2(2):end,1:end-d(3),:),struct('type','()','subs',{{':'}})), ...
                    subsref(X(1 + d2(1):end,1 + d2(2):end,1:end-d(3),:),struct('type','()','subs',{{':'}})), ...
                    subsref(X(1:end-d2(1),1:end-d2(2),1 + d(3):end,:),struct('type','()','subs',{{':'}})), ...
                    subsref(X(1 + d(1):end - d(1),1:end-d2(2),1 + d(3):end,:),struct('type','()','subs',{{':'}})), ...
                    subsref(X(1 + d2(1):end,1:end-d2(2),1 + d(3):end,:),struct('type','()','subs',{{':'}})), ...
                    subsref(X(1:end-d2(1),1 + d(2):end - d(2),1 + d(3):end,:),struct('type','()','subs',{{':'}})), ...
                    subsref(X(1 + d(1):end - d(1),1 + d(2):end - d(2),1 + d(3):end,:),struct('type','()','subs',{{':'}})), ...
                    subsref(X(1 + d2(1):end,1 + d(2):end - d(2),1 + d(3):end,:),struct('type','()','subs',{{':'}})), ...
                    subsref(X(1:end - d2(1),1 + d2(2):end,1 + d(3):end,:),struct('type','()','subs',{{':'}})), ...
                    subsref(X(1 + d(1):end - d(1),1 + d2(2):end,1 + d(3):end,:),struct('type','()','subs',{{':'}})), ...
                    subsref(X(1 + d2(1):end,1 + d2(2):end,1 + d(3):end,:),struct('type','()','subs',{{':'}}))],[],2);
            elseif all(layer.mpSize == [3, 3, 1])
                d2 = 2.*d;
                [Y,ind] = max([subsref(X(1:end-d2(1),1:end-d2(2),:,:),struct('type','()','subs',{{':'}})), ...
                    subsref(X(1 + d(1):end - d(1),1:end-d2(2),:,:),struct('type','()','subs',{{':'}})), ...
                    subsref(X(1 + d2(1):end,1:end-d2(2),:,:),struct('type','()','subs',{{':'}})), ...
                    subsref(X(1:end-d2(1),1 + d(2):end - d(2),:,:),struct('type','()','subs',{{':'}})), ...
                    subsref(X(1 + d(1):end - d(1),1 + d(2):end - d(2),:,:),struct('type','()','subs',{{':'}})), ...
                    subsref(X(1 + d2(1):end,1 + d(2):end - d(2),:,:),struct('type','()','subs',{{':'}})), ...
                    subsref(X(1:end - d2(1),1 + d2(2):end,:,:),struct('type','()','subs',{{':'}})), ...
                    subsref(X(1 + d(1):end - d(1),1 + d2(2):end,:,:),struct('type','()','subs',{{':'}})), ...
                    subsref(X(1 + d2(1):end,1 + d2(2):end,:,:),struct('type','()','subs',{{':'}}))],[],2);
            else
                %might be much slower and more memory intensive than the hard
                %coded pooling windows above. would be good to find a better
                %implementation here.
                [Y, ind] = layer.generalPooling(X);
            end
            sX = size(X);
            sX(end + 1:4) = 1;
            targetSize = layer.targetSize(sX);
            Y = reshape(Y, targetSize);
            layer.mpInd = reshape(ind,targetSize);
        end
        
        function [DEDX, DEDW] = backward(layer, ~, DEDY)
            d = layer.sparsityFactor;
            %size difference between input and output
            sD = d.*(layer.mpSize - 1);
            sDEDY = size(DEDY);
            sDEDY(end+1:4) = 1;
            DEDX = zeros([sDEDY(1:3) + sD,sDEDY(4)],'like',DEDY);

            %get global indices in DEDX for values in DEDY
            indConv = layer.indConversion(size(DEDX));
            zeroMat = reshape(1:numel(DEDX),size(DEDX));
            gl_ind = indConv(layer.mpInd) + zeroMat(1:end-sD(1),1:end-sD(2),1:end-sD(3),:);

            %determine which values need to be propagated back to the
            %same index in DEDX
            [u_gl_ind,~,ic] = unique(gl_ind);

            %sum values in DEDY which came from same position in DEDX
            DEDY_summed = accumarray(ic,DEDY(:));

            %assign summed DEDY to corresponding positions in DEDX
            DEDX(u_gl_ind) = DEDY_summed;

            DEDW = [];
        end
        
        function indConv = indConversion(layer, szY)
            d = layer.sparsityFactor;
            indConv = zeros(layer.mpSize);
            v = permute(0:layer.mpSize(1) - 1,[2, 1]);
            indConv = bsxfun(@plus,indConv,v.*d(1));
            v = 0:layer.mpSize(2) - 1;
            indConv = bsxfun(@plus,indConv,v.*d(2)*szY(1));
            v = permute(0:layer.mpSize(3) - 1, [1 3 2]);
            indConv = bsxfun(@plus,indConv,v.*d(3)*szY(1)*szY(2));
            indConv = indConv(:);
        end
        
        function [Y, ind] = generalPooling(layer, X)
            d = layer.sparsityFactor;
            szX = size(X);
            szX(end + 1:4) = 1;
            mpRegionIndices = layer.getPoolRegionIndices(szX);
            [Y, ind] = max(X(mpRegionIndices),[],2);
            Y = reshape(Y,[floor((szX(1:3) - d.*(layer.mpSize - 1))./layer.stride),szX(4)]);
            layer.mpInd = reshape(ind,[floor((szX(1:3) - d.*(layer.mpSize - 1))./layer.stride),szX(4)]);
        end
        
        function ind = getPoolRegionIndices(layer, szX)
            %add this to generalPooling function and hope it can be
            %optimized
            d = layer.sparsityFactor;
            mpSiz = layer.mpSize;
            dFilterSize = [szX(1:3) - d.*(mpSiz - 1),szX(4)];
            indStart = layer.indConversion(szX) + 1;
            ind = bsxfun(@plus,indStart,0:layer.stride(1):dFilterSize(1) - 1);
            ind = ind(:);
            ind = bsxfun(@plus,ind,0:layer.stride(2)*szX(1):szX(1)*(dFilterSize(2) - 1));
            ind = ind(:);
            ind = bsxfun(@plus,ind,0:layer.stride(3)*szX(1)*szX(2):szX(1)*szX(2)*(dFilterSize(3) - 1));
            ind = ind(:);
            ind = bsxfun(@plus,ind,0:szX(1)*szX(2)*szX(3):szX(1)*szX(2)*szX(3)*(szX(4) - 1));
            ind = reshape(ind,prod(mpSiz),numel(ind)/prod(mpSiz))';
        end
        
        function targetSiz = targetSize(layer, inputSize)
            targetSiz = [inputSize(1:3) - layer.sparsityFactor.* ...
                        (layer.mpSize - 1),inputSize(4)];
        end
    end
    
end
