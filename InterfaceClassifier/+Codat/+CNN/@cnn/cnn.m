classdef cnn
    %CNN implementation using only convolutions with strides. It uses
    %sparse convolutions to do pixel classification.
    %
    % Constructor arguments
    %   cLayer: Total number of convolutional layer (mlp-layer are through
    %           of as convolutions with a kernel of size 1) including the
    %           output layer (i.e. hidden layer and output layer).
    %   featureMaps: Array of length cLayer + 1 specifying the number of 
    %           featureMaps in each layer including the input and output
    %           layer (number of neurons for mlp layer). Feature maps in
    %           input and output layer correspond to channels.
    %   filterSize: Cell array containing size of convolutional filter in
    %               each hidden (cLayer) and the output layer in all three
    %               dimensions in an array. If only one filter size is
    %               specified than all layers use the same filter size.
    %               (e.g. filterSize = {[3, 3, 3],[1, 1, 1]} for cLayer = 2
    %                     or filterSize = {[3, 3, 3]} for an arbitrary
    %                     number of layers).
    %   stride: Cell array containing the strides of the convolution
    %           (i.e. distance between receptive filds of neighboring
    %           neurons) in each dimension (specified as filterSize). 
    %           Strides are applied to the output of the current layer.
    %           The stride cell array can either have the size cLayer, one
    %           dimension more to stride the input layer or one stride for
    %           all layers.
    %   dropout: Array specifying the dropout for each layer. Dropout is
    %           applied to the output of a layer. The size of the dropout
    %           array can be cLayer of cLayer + 1 if dropout should also be
    %           applied to the input layer. Dropout in the last layer will
    %           automatically be set to zero.
    %   maxPool: Binary array specifying whether to perform 2x2x2 pooling
    %            after convolution and non-linearity in each layer.
    %   shortcut: Array of length (cLayer + 1) which specifies defining a
    %             shortcut connection for each layer from the specified
    %             layer. 0 means no shortcut connection for the
    %             corresponding layer. E.g. [0 0 0 2] would mean that layer
    %             4 receives a shortcut connection from layer 2. Currently
    %             only one shortcut connection is allowed per layer.
    %             Furthermore, layers that are connected via a shortcut
    %             connection must have the same number of feature maps.
    %             Reference: He, Kaiming, et al. "Deep Residual Learning
    %               for Image Recognition." arXiv preprint arXiv:1512.03385
    %               (2015).
    %   batchNorm: Array of length cLayer of logicals specifying for
    %              each layer whether batch normalization should be
    %              performed. Batch norm for the first and last layer is
    %              enforced to be false.
    %              Reference: Ioffe, Sergey, and Christian Szegedy. "Batch
    %               normalization: Accelerating deep network training by
    %               reducing internal covariate shift." arXiv preprint
    %               arXiv:1502.03167 (2015).
    %   optimizer: A Codat.Optimizer object.
    %
    % Properties
    %   layer: Total number of layers (= cLayer + 1).
    %   featureMaps: see above (constructor)
    %   filterSize: see above (constructor)
    %   dFilterSize: Filter size after d-regularity.
    %   stride: see above (constructor)
    %   W: Cell array of weights in each layer of the form
    %      W{layer}(x,y,z,featureMaps(layer-1),featureMaps(layer))
    %   b: Cell array of biases. One bias per feature map.
    %   d: stacked d-regularity factor per layer.
    %   Wmask: Mask for d-regularity used in backprop.
    %   dropout: see above(constructor).
    %   shortcut: see above (constructor)
    %   batchNorm: see above (constructor)
    %   optimizer: Instance of Codat.Optimizer used to train cnet.
    %   l2WeightDecayLambda: L2 regularizer factor.
    %   nonLinearity: Cell array of function handles containing the
    %       non-linearity for each layer.
    %   nonLinearityD: Cell array of function handles containing the
    %       derivative of the corresponding non-linearity.
    %   lossFunction: String containing the loss function used for
    %       training. Possible loss functions are
    %           'Squared'
    %           'cross-entroy' - only with sigmoid in last layer
    %           'softmax' - only with softmax in last layer
    %   actvtClass: Anonymous function converting to desired precision
    %       (i.e. @single or @(x)gpuArray(single(x)) see also
    %       setParamsToActvtClass
    %   border: Border around each target pixel for cnn architecture.
    %   numParams: Total number of cnn parameters.
    %   isTraining: Flag used during training to apply dropout.
    %   convAlg: String specifying the way of calculation the convolutions.
    %           Options are
    %           'fft1' - Fastest but very memory consuming. Only
    %                    applicable for small target sizes. (usually ~2x
    %                    faster then fft2).
    %           'fft2' - FFT convolution. (usually the best to use for
    %                    prediction). Note that backpropagation is
    %                    currently only implemented for a mode similar to
    %                    fft1.
    %           'convn' - MATLAB convn convolution. (Very slow for
    %                     backprop)
    %   memoryLimit: Very rough limitation of memory consumption for mlp
    %           layers. (only to prevent huge memory consumption due to
    %           parallelization - DONT RELY ON THIS TO SOLVE ALL MEMORY
    %           RELATED PROBLEMS!)
    % Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>
    
    properties
        layer                   %total number of layer
        featureMaps             %feature maps per layer
        filterSize              %filter size per layer
        dFilterSize             %filter size after d-regularity
        stride                  %stride after each convolution
        maxPool                 %perform maxPool after non-linearity
        W                       %cnn weights
        b                       %cnn biases
        d                       %total d-regularity after layer (i.e.
                                %stride and max pooling)
        Wmask                   %weight masks for d-regulary
        dropout                 %percentage of points to drop out
        optimizer               %optimizer for parameter udpate
        l2WeightDecayLambda = 0 %l2 weight decay factor
        nonLinearity
        nonLinearityD
        lossFunction
        actvtClass = @single    %precision/gpuArray of quantities
        border                  %total network border
        numParams               %total number of network parameters
        isTraining              %used for dropout
        convAlg                 %algorithm used for conv in forward pass
        memoryLimit             %approximate RAM limit for some large array operatios in GB
        shortcut                %specify for each layer from which previous
                                %layer it receives a shortcut connection
        batchNorm               %application of batch normalization (bn)
        bn_beta
        bn_gamma
        bn_muInf
        bn_sig2Inf
        
    end
    
    methods
        %constructor
        function cnet = cnn(cLayer, featureMaps, filterSize, stride, maxPool, dropout, shortcut, batchNorm, nonLinearity, loss, optimizer)
            if length(featureMaps) ~= cLayer + 1
                error('Specify the number of feature maps for each layer.');
            end
            cnet.layer = cLayer + 1;
            cnet.featureMaps = featureMaps;
            
            if length(filterSize) == 1
                cnet.filterSize = [{[]},repmat(filterSize,1,cLayer)];
            elseif length(filterSize) == cLayer
                cnet.filterSize = [{[]},filterSize];
            else
                error('Wrong filter size specified');
            end
            
            if length(stride) == 1
                cnet.stride = [{[1 1 1]},repmat(stride,1,cLayer)];
            elseif length(stride) == cLayer
                cnet.stride = [{[1 1 1]},stride];
            elseif length(stride) == cLayer + 1
                cnet.stride = stride;
            else
                error('Wrong stride specified');
            end
            
            if length(dropout) == 1
                cnet.dropout = [repmat(dropout,1,cLayer) 0];
            elseif length(dropout) == cLayer + 1
                cnet.dropout = dropout;
                cnet.dropout(cLayer + 1) = 0;
            else
                error('Wrong dropout fraction specified');
            end
            
            if length(shortcut) == 1 && shortcut == 0
                cnet.shortcut = [0 repmat(shortcut,1,cLayer)];
            elseif length(shortcut) == cLayer + 1
                cnet.shortcut = shortcut;
            else
                error('Wrong shortcut specified');
            end
            
            if length(batchNorm) == 1
                cnet.batchNorm = [repmat(batchNorm,1,cLayer) false];
            elseif length(batchNorm) == cLayer
                cnet.batchNorm = [batchNorm, false];
            elseif length(shortcut) == cLayer + 1
                cnet.batchNorm = batchNorm;
                cnet.batchNorm(cLayer + 1) = false;
            else
                error('Wrong batchNorm specified');
            end
            cnet.batchNorm(1) = false;
            
            if length(maxPool) == 1
                cnet.maxPool = [false,repmat(maxPool,1,cLayer)];
            elseif length(maxPool) == cLayer
                cnet.maxPool = [false,maxPool];
            elseif length(maxPool) == cLayer + 1
                cnet.maxPool = maxPool;
            else
                error('Wrong max pool specified');
            end
            if cnet.maxPool(cnet.layer)
                error('The last layer is not allowed to to max-pooling.');
            end
            
            %set d-sparsity factor from stride and max pooling
            tmp = ones(length(cnet.maxPool),3);
            tmp(cnet.maxPool,:) = repmat([2,2,2],sum(cnet.maxPool),1);
            cnet.d = cumprod(cell2mat(cnet.stride').*tmp);
            
            %set non-linearity and loss
            cnet.nonLinearity = cell(cnet.layer,1);
            cnet.nonLinearityD = cell(cnet.layer,1);
            if ischar(nonLinearity)
                nonLinearity = repmat({nonLinearity},1,cnet.layer);
            elseif length(nonLinearity) == cLayer
                nonLinearity = [{[]},nonLinearity];
            end
            for lyr = 2:cnet.layer
                switch nonLinearity{lyr}
                    case 'linear'
                        cnet.nonLinearity{lyr} = @(x)x;
                        cnet.nonLinearityD{lyr} = @(x)1;
                    case 'tanh'
                        cnet.nonLinearity{lyr} = @Codat.CNN.cnn.tanh;
                        cnet.nonLinearityD{lyr} = @Codat.CNN.cnn.tanhD;
                    case 'sigmoid'
                        cnet.nonLinearity{lyr} = @Codat.CNN.cnn.sigmoid;
                        cnet.nonLinearityD{lyr} = @Codat.CNN.cnn.sigmoidD;
                    case 'relu'
                        cnet.nonLinearity{lyr} = @Codat.CNN.cnn.relu;
                        cnet.nonLinearityD{lyr} = @Codat.CNN.cnn.reluD;
                    case 'elu'
                        cnet.nonLinearity{lyr} = @Codat.CNN.cnn.elu;
                        cnet.nonLinearityD{lyr} = @Codat.CNN.cnn.eluD;
                    case 'softmax'
                        if lyr == cnet.layer
                            cnet.nonLinearity{lyr} = @Codat.CNN.cnn.softmax;
                        else
                            error('Softmax non-linearity only implemented in output layer');
                        end
                        if ~strcmp(loss,'softmax')
                            error('Softmax non-linearity required in last layer for softmax error function');
                        end
                    otherwise
                        error('Non-linearity not defined');
                end
            end
            if strcmp(loss,'cross-entropy')
                if ~strcmp(nonLinearity{end},'sigmoid')
                    error('Cross-entropy function only supported for sigmoid activation in last layer');
                end
            end
            cnet.lossFunction = loss;
            
            %initalize weights and bias
            for lyr = 2:cnet.layer
                factor = 1;
                if cnet.dropout(lyr - 1) > 0
                    factor = factor/sqrt(1 - cnet.dropout(lyr));
                end
                if any(strcmp(nonLinearity{lyr},{'tanh','sigmoid','softmax'}))
                    factor = factor/sqrt(prod(cnet.filterSize{lyr}))/sqrt(cnet.featureMaps(lyr - 1));
                elseif any(strcmp(nonLinearity{lyr},{'relu','elu'}))
                    factor = factor*sqrt(2)/sqrt(prod(cnet.filterSize{lyr}))/sqrt(cnet.featureMaps(lyr - 1));
                end
                if cnet.shortcut(lyr) > 0
                    factor = factor/sqrt(2);
                end
                tmp = factor.*randn([cnet.filterSize{lyr}, cnet.featureMaps(lyr - 1), cnet.featureMaps(lyr)],'single');
                cnet.W{lyr} = tmp;
                [~, cnet.Wmask{lyr}] = cnet.sparseKernel(tmp,cnet.d(lyr - 1,:));
                cnet.dFilterSize{lyr} = cnet.mSize(cnet.Wmask{lyr},1:3);
                cnet.b{lyr} = zeros(cnet.featureMaps(lyr),1,'single');
            end
            
            %initialize bn parameters
            for lyr = 2:cnet.layer
                if cnet.batchNorm(lyr)
                    cnet.bn_beta{lyr} = zeros(1,1,1,cnet.featureMaps(lyr));
                    cnet.bn_gamma{lyr} = ones(1,1,1,cnet.featureMaps(lyr));
                    cnet.bn_muInf{lyr} = zeros(1,1,1,cnet.featureMaps(lyr));
                    cnet.bn_sig2Inf{lyr} = ones(1,1,1,cnet.featureMaps(lyr));
                end
            end
            
            cnet.border = cnet.calculateBorder();
            cnet.numParams = cnet.calculateNumParams();
            cnet.l2WeightDecayLambda = 0;
            cnet.isTraining = false;
            cnet.convAlg = 'fft1';
            cnet.memoryLimit = 6;
            cnet.optimizer = optimizer.init(cnet.numParams);
        end
        
        function border = calculateBorder(cnet)
            border = 0;
            for l = 2:cnet.layer
                % border coming from stride and max pooling
                border = border + cnet.d(l - 1,:).*(cnet.filterSize{l} - 1) + cnet.maxPool(l).*cnet.stride{l}.*cnet.d(l - 1,:);
            end
            if floor(border/2) ~= border/2
                error('Border of [%s] is not symmetric around target points.',num2str(border));
            end
        end
        
        function vec = param2Vec(cnet,varargin)
            if isempty(varargin)
                PStruct.W = cnet.W;
                PStruct.b = cnet.b;
                PStruct.beta = cnet.bn_beta;
                PStruct.gamma = cnet.bn_gamma;
            else
                PStruct = varargin{1};
            end
            vec = zeros(0,'like',cnet.W{2});
            switch cnet.convAlg
                case {'fft1','fft2'}
                    for l = 2:cnet.layer
                        vec = cat(1,vec,PStruct.W{l}(:));
                        vec = cat(1,vec,PStruct.b{l});
                        if cnet.batchNorm(l)
                            vec = cat(1,vec,PStruct.beta{l}(:),PStruct.gamma{l}(:));
                        end
                    end
                case 'convn'
                    for l = 2:cnet.layer
                        vec = cat(1,vec,shiftdim(PStruct.W{l}(repmat(cnet.Wmask{l},[1, 1, 1, cnet.featureMaps(l - 1), cnet.featureMaps(l)]))));
                        vec = cat(1,vec,PStruct.b{l});
                        if cnet.batchNorm(l)
                            vec = cat(1,vec,PStruct.beta{l}(:),PStruct.gamma{l}(:));
                        end
                    end
            end
        end
        
        function cnet = vec2Param(cnet, vec)
            count = 1;
            switch cnet.convAlg
                case {'fft1','fft2'}
                    for l = 2:cnet.layer
                        %fill weights
                        nl = [cnet.filterSize{l}, cnet.featureMaps(l - 1), cnet.featureMaps(l)];
                        tmpCount = prod(nl);
                        cnet.W{l} = reshape(vec(count: count + tmpCount - 1),nl);
                        count = count + tmpCount;
                        %fill bias
                        tmpCount = length(cnet.b{l});
                        cnet.b{l} = vec(count:count + tmpCount - 1);
                        count = count + tmpCount;
                        %fill bn
                        if cnet.batchNorm(l)
                            cnet.bn_beta{l} = reshape(vec(count:count + cnet.featureMaps(l) - 1),1,1,1,cnet.featureMaps(l));
                            count = count + cnet.featureMaps(l);
                            cnet.bn_gamma{l} = reshape(vec(count:count + cnet.featureMaps(l) - 1),1,1,1,cnet.featureMaps(l));
                            count = count + cnet.featureMaps(l);
                        end
                    end
                case 'convn'
                    for l = 2:cnet.layer
                        %fill weights
                        nl = [cnet.filterSize{l}, cnet.featureMaps(l - 1), cnet.featureMaps(l)];
                        tmpCount = prod(nl);
                        cnet.W{l}(repmat(cnet.Wmask{l},[1, 1, 1, cnet.featureMaps(l - 1), cnet.featureMaps(l)])) = reshape(vec(count: count + tmpCount - 1),nl);
                        count = count + tmpCount;
                        %fill bias
                        tmpCount = length(cnet.b{l});
                        cnet.b{l} = vec(count:count + tmpCount - 1);
                        count = count + tmpCount;
                        %fill bn
                        if cnet.batchNorm(l)
                            cnet.bn_beta{l} = reshape(vec(count:count + cnet.featureMaps(l) - 1),1,1,1,cnet.featureMaps(l));
                            count = count + cnet.featureMaps(l);
                            cnet.bn_gamma{l} = reshape(vec(count:count + cnet.featureMaps(l) - 1),1,1,1,cnet.featureMaps(l));
                            count = count + cnet.featureMaps(l);
                        end
                    end
            end
        end
        
        function numParams = calculateNumParams(cnet)
            %weights
            numParams = sum(cell2mat(cellfun(@(x)prod(x),cnet.filterSize(2:end),'UniformOutput',false)).*cnet.featureMaps(2:end).*cnet.featureMaps(1:end-1));
            %biases
            numParams = numParams + sum(cellfun(@(x)numel(x),cnet.b));
            %bn params
            numParams = numParams + 2*sum(cnet.batchNorm.*cnet.featureMaps);
        end
        
        function cnet = saveobj(cnet)
            %convert function handles to string
            cnet.nonLinearity{1} = [];
            cnet.nonLinearityD{1} = [];
            cnet.nonLinearity(2:end) = cellfun(@func2str,cnet.nonLinearity(2:end),'UniformOutput',false);
            cnet.nonLinearityD(2:end) = cellfun(@func2str,cnet.nonLinearityD(2:end),'UniformOutput',false);
        end
    end
    
    methods (Static)
        function Y = flipdims(X)
            Y = reshape(X(end:-1:1),size(X));
        end
        
        [Wd,mask] = sparseKernel(W,d)
        
        function X = tanh(X)
            X = 1.7159*tanh(0.66*X);
        end
        
        function X = tanhD(X)
            X = 0.66*(1.7159 - 1./1.7159.*X.^2);
        end
        
        function X = sigmoid(X)
            X = 1./(1+exp(-X));
        end
        
        function X = sigmoidD(X)
            X = X.*(1 - X);
        end
        
        function X = relu(X)
            X = cast(max(0,X),'like',X);
        end
        
        function X = reluD(X)
            X = cast(X > 0,'like',X);
        end
        
        function X = elu(X)
            X(X < 0) = exp(X(X < 0)) - 1;
        end
        
        function X = eluD(X)
            X = cast(X >= 0,'like',X) + (X + 1).*(X < 0);
        end
        
        y = softmax(x);
        
        function d = mSize(A,dims)
            %return sizes for specified dims only
            d = size(A);
            d(end + 1:max(dims)) = 1;
            d = d(dims);
        end
        
        y = fftd( A, siz, dim, d);
        
        A = cropActivation(A,targetSize);
        A = padArray(A,targetSize);
        
        function cnet = loadobj(cnet)
            %convert non-linearity strings back to function handles
            cnet.nonLinearity(2:end) = cellfun(@str2func,cnet.nonLinearity(2:end),'UniformOutput',false);
            cnet.nonLinearityD(2:end) = cellfun(@str2func,cnet.nonLinearityD(2:end),'UniformOutput',false);
        end
    end
    
end

