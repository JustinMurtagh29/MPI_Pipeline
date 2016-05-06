classdef denseConvolutionalLayer < Codat.Net.Layer.baseLayer
    %DENSECONVOLUTIONALLAYER 3D convolutional layer for dense pixel
    % classification.
    %
    % PROPERTIES
    % filterSize: Integer vector of length 3 specifying the filter size in
    %             each dimension.
    % featureMapsOut: The number of output feature maps/output channels.
    % featureMapsIn: The number of input feature maps/input channels.
    % sparsityFactor: Integer vector of length 3 defining the sparsity of the
    %         pooling window. This is only important for dense pixel
    %         classification and depends on the strides of previous layers.
    %         The sparsity factor of this layer should be the product of
    %         stides of preceding layers.
    % W: Array of double containing the layer weight parameters. The size
    %    of W is [filterSize, featureMapsIn, featureMapsOut].
    % b: Vector of double containing the layer bias parameters. b contains
    %    one bias for each output feature map.
    % convAlg: String specifying which implementation to use for the
    %          convolution. All implementations produce the same output but
    %          differ in speed and memory requirement.
    %          Options are
    %          'fft': Convolution via fft and full vectorization.
    %                     This is typically faster but more memory
    %                     consuming than 'convn'.
    %          'fft2': Another fft implementation of convolutions requiring
    %                  less memory. This mode uses the same routine as 'fft' for
    %                  backpropagation.
    %          'convn': Convolution via MATLABs convn function. Note that this
    %                   mode was not optimized.
    % propDown: Flag indicating whether the error w.r.t. the input should be
    %           calculated during backpropagation.
    % memoryLimit: Double trying to keep the memory consumption for 1x1x1
    %              convolutions below that value (in GB).
    % wMask: Mask used to deal with sparse filter for convAlg 'fft2' and
    %        'convn'.
    %
    % Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

    properties
        filterSize
        sparseFilterSize
        featureMapsOut
        featureMapsIn
        sparsityFactor
        W
        b
        convAlg = 'fft';
        propDown = true;
        memoryLimit = 8;
        wMask;
    end

    methods
        function obj = denseConvolutionalLayer(filterSize, stride, featureMapsIn, featureMapsOut, sparsityFactor)
            % Constructor. See clas properties for description of inputs.

            obj.filterSize = filterSize;
            obj.stride = stride;
            obj.featureMapsOut = featureMapsOut;
            obj.featureMapsIn = featureMapsIn;
            obj.numParams = prod(filterSize)*featureMapsIn*featureMapsOut + featureMapsOut;
            if ~exist('sparsityFactor','var') || isempty(sparsityFactor)
                obj.sparsityFactor = [1, 1, 1];
            else
                obj.sparsityFactor = sparsityFactor;
            end
            obj.sparseFilterSize = filterSize + (filterSize - 1).*(sparsityFactor - 1);
            obj.W = zeros([filterSize, featureMapsIn, featureMapsOut],'single');
            obj.b = zeros(obj.featureMapsOut,1,'single');
            obj.wMask = true(obj.filterSize);
        end

        function [layer, Y] = forward(layer, X)
            % Forward pass through layer.
            % INPUT X: 4d array where the first three dimensions correspond to
            %          spatial coordinates and the fourth to channels/feature
            %          maps.

            sX = size(X);
            sX(end+1:4) = 1;
            Y = zeros([sX(1:3) - layer.sparseFilterSize + 1, layer.featureMapsOut], 'like', X);
            %mlpConv layer (do direct multiplication instead of fft/convn)
            if prod(layer.filterSize) == 1  && numel(X)*layer.featureMapsOut*4 < layer.memoryLimit*1073741824
                Y = bsxfun(@plus,squeeze(sum(bsxfun(@times,X,flip(layer.W,4)),4)),permute(layer.b,[4 3 2 1]));
            elseif prod(layer.filterSize) == 1
                for fm = 1:layer.featureMapsOut
                    Y(:,:,:,fm) = sum(bsxfun(@times,X,flip(layer.W(:,:,:,:,fm),4)),4) + layer.b(fm);
                end
            else %convolutional layer
                switch layer.convAlg
                    case 'fft' %semi vectorized and adapted for strided convolutions
                        %valid convolution via fft over whole array (high memory
                        %requirement)

                        w = flip(layer.W,4);
                        %determine size for fft (should be multiple of
                        %sparsityFactor in all dimensions)
                        d = layer.sparsityFactor;
                        ffSize = sX(1:3) + mod(-sX(1:3),d);
                        wSize = layer.sparseFilterSize;

                        %only do fft if filter size > 1
                        fftDims = layer.filterSize > 1;
                        ffSize(~fftDims) = sX(~fftDims);
                        fftDims = find(fftDims);
                        for k = fftDims
                            X = fft(X,ffSize(k),k);
                            w = layer.fftd(w,ffSize(k),k, d(k));
                        end

                        %fully vectorized product in fourier space
                        %             C = bsxfun(@times,activity,W);
                        %             C = permute(sum(C,4),[1, 2, 3, 5, 4]);

                        %semi-vectorized with usually much less memory requirement and
                        %almost same speed
                        C = zeros([ffSize(1:3), layer.featureMapsOut],'like',X);
                        for i = 1:layer.featureMapsOut
                            C(:,:,:,i) = sum(bsxfun(@times,X,w(:,:,:,:,i)),4);
                        end

                        for k = fftDims
                            C = ifft(C,ffSize(k),k);
                        end
                        %make sure is real and get indices for valid conv
                        Y = bsxfun(@plus,real(C(wSize(1):sX(1), ...
                            wSize(2):sX(2),wSize(3):sX(3),:)), ...
                            permute(layer.b,[4 3 2 1]));

%                    case 'fft2'
%                        w = flip(layer.W,4);
%                        sW = layer.sparseFilterSize;
%                        fftDims = layer.filterSize > 1;
%                        d = layer.sparsityFactor;
%                        ffSize = sX(1:3) + mod(-sX(1:3),d);
%                        wSize = layer.sparseFilterSize;
%                        for k = fftDims
%                            X = X = fft(X,ffSize(k),k);
%                        end
%                        for fm = 1:layer.featureMapsOut
%                            for k = fftDims
%                                wFm = layer.fftd(w(:,:,:,:,fm),ffSize(k),k, d(k));
%                            end
%                            Y(:,:,:,fm) = subsref(real(ifftn(X.*wFM)), ...
%                                struct('type','()','subs', ...
%                                {{sW(1):sX(1),sW(2):sX(2),sW(3):sX(3), ...
%                                layer.featureMapsIn}})) + layer.b(fm);
%                        end

                    case 'fft2'
                        w = layer.W;
                        sW = layer.sparseFilterSize;
                        X = fftn(X);
                        for fm = 1:layer.featureMapsOut
                            wFM = fftn(w(:,:,:,:,fm),sX);
                            Y(:,:,:,fm) = subsref(real(ifftn(X.*wFM)), ...
                                struct('type','()','subs', ...
                                {{sW(1):sX(1),sW(2):sX(2),sW(3):sX(3), ...
                                layer.featureMapsIn}})) + layer.b(fm);
                        end

                    case 'convn'
                        for fm = 1:layer.featureMapsOut
                            Y(:,:,:,fm) = convn(X,layer.W(:,:,:,:,fm),'valid') ...
                                + layer.b(fm);
                        end
                end
            end
        end

        function [DEDX, DEDW] = backward(layer, X, DEDY)
            %Error backpropagation through layer.
            % INPUT X: Input array from forward pass.
            %       DEDY: Derivative of error function wrt to output Y of layer.

            DEDX = zeros(size(X),'like',X);
            if prod(layer.filterSize) == 1 %mlp layer
                %bias gradient
                gb = squeeze(sum(sum(sum(DEDY,1),2),3));
                %weight gradient
                if numel(X)*size(DEDY,4)*4 < layer.memoryLimit*1073741824
                    gW = sum(sum(sum(bsxfun(@times,flip(X,4),permute(DEDY,[1 2 3 5 4])),1),2),3);
                    if layer.propDown
                        DEDX = sum(bsxfun(@times,permute(DEDY,[1 2 3 5 4]),flip(layer.W{lyr},4)),5);
                    end
                else
                    gW = zeros(size(layer.W),'like',layer.W);
                    for fm = 1:layer.featureMaps(lyr)
                        gW(:,:,:,:,fm) = layer.flipdims(sum(sum(sum(bsxfun(@times,X,DEDY(:,:,:,fm)),1),2),3));
                        if layer.propDown
                            DEDX = DEDX + bsxfun(@times,DEDY(:,:,:,fm),layer.flipdims(layer.W{lyr}(:,:,:,:,fm)));
                        end
                    end
                end
            else %convolutional layer
                gb = shiftdim(sum(sum(sum(DEDY,1),2),3));
                switch layer.convAlg
                    case 'fft'
                        %transform to fourier space
                        actSize = size(X);
                        actSize(end+1:4) = 1;
                        d = layer.sparsityFactor;
                        ffSize = actSize(1:3) + mod(-actSize(1:3),d);
                        delSize = size(DEDY);
                        delSizeNew = [ffSize(1:3), size(DEDY,4)];
                        fftDimIdx = layer.filterSize > 1;
                        ffSize(~fftDimIdx) = actSize(~fftDimIdx);
                        fftDims = find(fftDimIdx);
                        X = layer.flipdims(X);
                        for dim = fftDims
                            X = fft(X,ffSize(dim),dim);
                            DEDY = fft(DEDY,delSizeNew(dim),dim);
                        end
                        for dim = setdiff(1:3,fftDims)
                            X = flip(X,dim);
                        end
                        gW = bsxfun(@times,X,permute(DEDY,[1 2 3 5 4]));

                        %gW calculation requires strided convolution with stride d
                        %which should be implemented by calculating ifft with stride d
                        %which is not natively implemented by matlab and hence full

                        %ifft is calculated and unnecessary indices discarded)
                        for dim = setdiff(1:3,fftDims) %accumulate gradient over batches
                            gW = sum(gW,dim);
                        end
                        fSize = layer.filterSize;

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
                        gW = real(gW);

                        if layer.propDown
                            w = flip(flip(flip(flip(layer.W,1),2),3),4);
                            DEDY = permute(DEDY,[1, 2, 3, 5, 4]);
                            for dim = fftDims
                                w = layer.fftd(w,ffSize(dim),dim,d(dim));
                            end
                            DEDX = sum(bsxfun(@times,DEDY,w),5);
                            for dim = fftDims
                                DEDX = ifft(DEDX,ffSize(dim),dim);
                            end
                            DEDX = real(DEDX(1:actSize(1),1:actSize(2),1:actSize(3),:));
                        end

                    case 'fft2'
                        gW = zeros([layer.sparseFilterSize, layer.featureMapsIn, layer.featureMapsOut],'like',layer.W);
                        for fm = 1:layer.featureMapsOut
                            gW(:,:,:,:,fm) = conv3fft(layer.flipdims(X),DEDY(:,:,:,fm),'valid');
                            if layer.propDown
                                DEDX = DEDX + conv3fft(DEDY(:,:,:,fm),layer.flipdims(layer.W(:,:,:,:,fm)));
                            end
                        end
                        %sparsify gradient (should be replaced by strided
                        %convolution)
                        gW = gW(repmat(layer.wMask,[1, 1, 1, layer.featureMapsIn, layer.featureMapsOut]));

                    case 'convn'
                        gW = zeros(size(layer.W),'like',layer.W);
                        for fm = 1:layer.featureMapsOut
                            gW(:,:,:,:,fm) = convn(layer.flipdims(X),DEDY(:,:,:,fm),'valid');
                            if layer.propDown
                                DEDX = DEDX + convn(DEDY(:,:,:,fm),layer.flipdims(layer.W(:,:,:,:,fm)));
                            end
                        end
                        %sparsify gradient (should be replaced by strided
                        %convolution)
                        gW = gW(repmat(layer.wMask,[1, 1, 1, layer.featureMapsIn, layer.featureMapsOut]));
                end
            end
            DEDW = [gW(:);gb(:)];
        end

        function layer = initializeWeights(layer, method, varargin)
            %Initialize the layer weights using the given method and
            %parameters.
            % INPUT method: String specifying the initialization method.
            %               Options are
            %               'fan-in': Initialize weights based on input size
            %                         to a single neuron.
            %               'randn': Initialization from gauss distribution.
            %       initParams: Parameters needed for weight initialization.
            %                   The needed parameters are
            %           'fan_in':
            %               1) The directly preceeding dropout fraction.
            %               2) The directly succeeding non-linearity as string.
            %           'randn': Standard deviation of the gaussian.

            switch method
                case 'fan-in'
                    factor = 1;
                    if varargin{1} > 0 && varargin{1} < 1
                        factor = factor/sqrt(1 - varargin{1});
                    end
                    if any(strcmp(varargin{2},{'tanh','sigmoid','softmax'}))
                        factor = factor/sqrt(prod(layer.filterSize))/sqrt(layer.featureMapsIn);
                    elseif strcmp(varargin{2},'relu')
                        factor = factor*sqrt(2)/sqrt(prod(layer.filterSize))/sqrt(layer.featureMapsIn);
                    end
                    tmp = factor.*randn([layer.filterSize, layer.featureMapsIn, layer.featureMapsOut],'single');
                    layer = layer.setWeights(tmp);
                case 'randn'
                    tmp = varargin{1}.*randn([layer.filterSize, layer.featureMapsIn, layer.featureMapsOut],'single');
                    layer = layer.setWeights(tmp);
                otherwise
                    error('Initialization method ''%s'' not implemented.',method);
            end
        end

        function layer = initializeBias(layer, method, varargin)
            %Initialize the layer bias using the given method and
            %parameters.
            % INPUT method: String specifying the initialization method.
            %               Options are
            %               'const': Initialize all bias parameters with
            %                   the same constant.
            %       varargin: Parameters needed for weight initialization.
            %                   The needed parameters are
            %           'const': The value for initialization.

            switch method
                case 'const'
                    layer.b(:) = varargin{1};
                otherwise
                    error('Initialization method ''%s'' not implemented.',method);
            end
        end

        function layer = castParams(layer, type)
            %Case laye rparameters to specified type.
            % INPUT type: Anonymous function which is applied to the layer
            %             parameters, e.g. @double, @gpuArray
            % NOTE This function always does gather before case params to avoid
            %      confusion when users try to do @double on a gpuArray.

            layer.W = type(gather(layer.W));
            layer.b = type(gather(layer.b));
        end

        function layer = setConvAlg(layer, mode)
            % INPUT mode: String specifying the convAlg.
            switch mode
                case 'fft'
                    if any(strcmp({'fft2','convn'},layer.convAlg))
                        layer.W = reshape(layer.W(repmat(layer.wMask,[1 1 1 layer.featureMapsIn layer.featureMapsOut])),[layer.filterSize layer.featureMapsIn layer.featureMapsOut]);
                        layer.wMask = true(layer.filterSize);
                    end
                case {'convn', 'fft2'}
                    if strcmp('fft',layer.convAlg)
                        [layer.W, layer.wMask] = layer.sparseKernel(layer.W,layer.sparsityFactor);
                    end
                otherwise
                    error('Convolutional algorithm ''%s'' not implemented',mode);
            end
            layer.convAlg = mode;
        end

        function params = param2Vec(layer)
            %Returns a vector containing all trainable layer parameters.

            if strcmp(layer.convAlg,'fft')
                params = [layer.W(:); layer.b];
            else
                weights = layer.getWeights();
                params = [weights(:); layer.b];
            end
        end

        function layer = vec2Param(layer, vec)
            % Set layer parameter to specified input vector.
            % INPUT vec: A vector of length(layer.numParams). Weights will be
            %            filled first along columns (i.e. like W(:) = vec) and
            %            then the bias.

            if length(vec) == layer.numParams
                wSize = [layer.filterSize, layer.featureMapsIn, layer.featureMapsOut];
                numWeights = prod(wSize);
                layer = layer.setWeights(vec(1:numWeights));
                layer.b = vec(numWeights + 1:end);
            else
                error('Input vector length does not match the number of parameters of this layer.');
            end
        end

        function W = getWeights(layer)
            % Returns a vector containing the layer weights.
            W = layer.W(repmat(layer.wMask,[1, 1, 1, layer.featureMapsIn, layer.featureMapsOut]));
        end

        function layer = setWeights(layer, W)
            %Set the layer weights to specified input vector.
            % INPUT W: Any array containing exactly the number of free weight
            %          parameters
            %          (i.e. prod(filtersize)*featureMapsIn*featureMapsOut).

            layer.W(repmat(layer.wMask,[1, 1, 1, layer.featureMapsIn, layer.featureMapsOut])) = W(:);
        end
    end

    methods (Static)
        function out = fftd( A, siz, dim, d )
            %FFTD Perform fft with d-regular array.
            % INPUT A: Input array.
            %       siz: Size of target in dim. Size is adapted to be processed
            %            effectively by fft so check on output size.
            %       dim: Apply fft along dim.
            %       d: D-regularity factor.
            % OUTPUT out: Fourier transformation of input array.

            siz = siz + mod(-siz,d);
            reSize = ones(1,ndims(A));
            reSize(dim) = d;
            tmp = fft(A,siz/d,dim);
            out = repmat(tmp,reSize);
        end

        function X = flipdims(X)
            X = reshape(X(end:-1:1),size(X));
        end

        function [Wd,mask] = sparseKernel(W, d)
            %SPARSEKERNEL Create d-regular sparse kernel from W.
            % INPUT W: 4D array, where the fourth dimension corresponds to
            %          feature maps and will not be transformed.
            %       d: Vector containing the regularity for the first three
            %          dimensions in W. (d = [2 2 2] means that every second
            %          row, colum, z-plane will be zeros.

            fSize = size(W);
            fdSize = fSize;
            if length(fSize) < 3
                fSize(3) = 1;
            end
            fSize = fSize(1:3);
            fsSize = fSize + (fSize - 1).*(d - 1);
            fdSize(1:3) = fsSize;
            Wd = zeros(fdSize,'like',W);
            xg = false(fsSize(1),1);
            yg = false(fsSize(2),1);
            zg = false(fsSize(3),1);
            xg(1:d(1):end) = true;
            yg(1:d(2):end) = true;
            zg(1:d(3):end) = true;
            [x,y,z] = meshgrid(xg,yg,zg);
            mask = x & y & z;
            repSize = fdSize;
            repSize(1:3) = 1;
            Wd(repmat(mask,repSize)) = W;
        end

        function C = conv3fft(A, B, shape)
            %CONV3FFT (Batch) 3d convn using fft.
            % Perform convolution using fft. This is much faster for large arrays.
            % This function is intended to be used for convolving two 3d array or two
            % 4d arrays where the fourth dimension has the same size ("batch 3d").
            % The fourth dimension should be a a power of 2.
            %
            % Note: The first time you execute this function MATLAB searches parameter
            % settings and saves the optimal ones internally which is rather slow.
            % Performance should be tested after this step was performed or previous
            % turning parameter should be loaded (see fftw).
            %
            % Recommended usage:
            %   3D: size(A),size(B) > 5.
            %
            % Warning: There are no sanity checks. Better run convn first to see
            % whether all your inputs are fine (i.e. sizes are ok....). Also A and B
            % need to have the same number of dimensions.
            %
            % INPUT A: N-d array.
            %       B: N-d array.
            %       shape: As for convn.

            if nargin < 3 || isempty(shape)
                shape = 'full';
            end

            switch shape
                case 'full'
                    ifun = @(m,n) 1:m+n-1;
                case 'valid'
                    ifun = @(m,n) n:m;
                case 'same'
                    ifun = @(m,n) ceil((n-1)/2)+(1:m);
                otherwise
                    error('conv3fft: unknown shape %s', shape);
            end

            %set fft size and indices for each dimension
            nd = max(ndims(A),ndims(B));
            subs(1:nd) = {':'};
            l = ones(1,4);
            for dim = 1:nd
                m = size(A,dim);
                n = size(B,dim);
                subs{dim} = ifun(m,n);
                %determine fft size
                switch shape
                    case 'full'
                        l(dim) = fftwn(m+n-1);
                    case 'valid'
                        l(dim) = fftwn(m);
                    case 'same'
                        l(dim) = fftwn(m + ceil((n-1)/2));
                end
            end

            A = fftn(A,l);
            B = fftn(B,l);
            C = ifftn(A.*B);
            C = real(C(subs{:}));
        end
    end

end
