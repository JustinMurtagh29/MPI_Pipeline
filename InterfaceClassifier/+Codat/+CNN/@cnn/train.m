function [ cnet, val_loss ] = train( cnet, stacksT, stacksV, options  )
%TRAIN Main CNN training function.
% Perform training on previously saved cubes. Cubes can futher be split
% into smaller cubes to reduce I/O overhead and augmented during training.
% StackT files and subcubes are randomly permuted during training.
%
% INPUT cnet: A Codat.CNN object.
%       stacksT: Cell array containing the paths to training files. Each
%                training file should contain the variables "raw" and
%                "target".
%       stacksV: Like stacksT for validation.
%       options: Training option struct containing
%           display:    Display information frequency (in iterations).
%           tr_size:    Size of target cube during training.
%           val_size:   Size of target cube during validation.
%           val_iter:   Validation frequency ( in iterations). (0 for no
%                       validation).
%           val_fwd_alg:Algorithm for forward pass during validation.
%           gpuDev:     Integer specifying the gpuDevice. (True corresponds
%                       to gpuDevice(1), false or 0 for no GPU training)
%           max_iter:   Maximum number of iterations. Specify Inf to sweep
%                       once through the whole training data. See also
%                       lr_policy "poly".
%           lr_policy: String specifying a learning rate policy
%               'constant': Learning rate does not change
%               'step': Adapt learning rate by factor (gamma) in every X
%                       iterations (step_size). Specify gamma and step_size
%                       for this lr_policy
%               'exp': Exponential decay of learning rate with the number
%                      of iterations by a gamma factor. Specify gamma for
%                      this lr_policy.
%               'poly': Polynomial decay of learning rate such that the
%                       learning is zero at max_iter. Specify gamma and
%                       max_iter for this lr_policy.
%           augment:    Bool whether to augment training set (factor 8x
%                       increase in data).
%           snapshot:   Snapshot frequency (in iterations). (0 for no
%                       snapshot)
%           snapshot_name: Path and name of snapshot file (without
%                       file-ending, iterations will be specified in
%                       addition). If empty than no snapshots will be
%                       saved.
%           save_imp: Bool whether to check for improvement at each
%                     validation and save a snapshot if current validation
%                     loss is the best.
%           data_pre: Function handle for data loading function as string.
%                     (use this for preprocessing, if empty the standard
%                     load will be used)
%           class_ratio: Maximal ratio between number of points of
%                     different classes. (Empty for no balancing).
%           skip_single_class_cubes: Bool indicating whether to skip cubes
%                      with only one class/target value.
%           gamma: See lr_policy "step" and "exp".
%           step_size: See lr_policy "step"
%           power: See lr_policy "poly".
%           rot_inv: Train on all 90 degree rotations of data to learn
%               rotation invariant features.
% OUTPUT cnet: Trained cnet.
%        val_loss: Validation loss during training.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

fprintf('[%s] Training CNN.\n',datestr(now));
%set options
actvtClass = cnet.actvtClass; %to cast back at end of training
if options.gpuDev
    fprintf('[%s] Training on gpu.\n',datestr(now));
    %gpuDevice(double(options.gpuDev)); %not needed anymore due to gpu exclusive mode 
    cnet = cnet.setParamsToActvtClass(@gpuArray);
else
    fprintf('[%s] Training on cpu.\n',datestr(now));
end
if isempty(options.data_pre)
    options.data_pre = @loadData;
    fprintf('[%s] Using standard load function.\n',datestr(now));
else
    fprintf('[%s] Using %s for data loading.\n',datestr(now), options.data_pre);
    options.data_pre = str2func(options.data_pre);
end
if ~options.augment
    options.augmentationIterations = 1;
    fprintf('[%s] Data augmentation set to false.\n',datestr(now));
else
    options.augmentationIterations = 8;
    fprintf('[%s] Using data augmentation during training.\n',datestr(now));
end
if ~isempty(options.class_ratio)
    fprintf('[%s] Class ratio balancing set to %d.\n',datestr(now),options.class_ratio);
end
if options.rot_inv
    fprintf('[%s] Using rotation invariant learning.\n',datestr(now));
end
fprintf('[%s] Learning rate policy set to %s.\n',datestr(now),options.lr_policy);

%training
rng('shuffle');
switch options.lr_policy
    case {'step','exp','poly'}
        base_lr = cnet.optimizer.learningRate;
end
total_iter = 1;
tr_loss = [];
val_loss = [];
skipVal = false;
options.trainStart = datestr(now);

%iterate over cubes
while total_iter <= options.max_iter
    for s = randperm(length(stacksT)) %all stack files
        fprintf('[%s] Training on stack file %s (%d/%d).\n',datestr(now),stacksT{s},s,length(stacksT));
        [input, target, targetWeights] = options.data_pre(stacksT{s});
        
        %perform augmentation if specified
        for augIt = 1:options.augmentationIterations
            [input, target, targetWeights] = Codat.CNN.Misc.invarianceOperation(augIt,input,target,targetWeights);
            iterations = prod(floor(cnet.mSize(target,1:3)./options.tr_size));
            rndIterIndices = randperm(iterations);
            
            if options.gpuDev
                %this is done here in order to minimize gpu transformation
                %overhead
                input = gpuArray(input);
                target = gpuArray(target);
                targetWeights = gpuArray(targetWeights);
            end
            
            %iterate over subcubes
            for iter = 1:iterations
                
                %validation
                if options.val_iter > 0 && mod(total_iter - 1, options.val_iter) == 0 && ~skipVal
                    fprintf('[%s] Iteration %d: Validating net.\n',datestr(now),total_iter - 1);
                    if options.rot_inv
                        val_loss(end + 1) = cnet.rotValidate(stacksV, options);
                    else %standard training
                        val_loss(end + 1) = cnet.validate(stacksV, options);
                    end
                    fprintf('[%s] Iteration %d: Validation loss: %.3f\n',datestr(now),total_iter - 1,val_loss(end));
                    if options.save_imp && total_iter > 1 && min(val_loss) == val_loss(end)
                        saveCnet(cnet, options, stacksT, stacksV, total_iter, val_loss);
                    end
                end
                
                %cube tiling according to train_size
                [x_train,y_train,curr_weights] = Codat.CNN.Misc.tileTrainingCubes(input, target, ...
                                    options.tr_size, cnet.border, options.max_iter, ...
                                    rndIterIndices(iter),targetWeights);
                if options.skip_single_class_cubes
                    if length(unique(y_train{1}(:))) == 1
                        skipVal = true;
                        continue
                    end
                end
                
                %class balancing
                if ~isempty(options.class_ratio)
                    curr_weights{1} = Codat.CNN.Misc.balanceClasses(y_train{1}, curr_weights{1}, options.class_ratio);
                    if all(curr_weights{1}(:))
                        skipVal = true;
                        continue
                    end
                end
                
                %actual cnn training iteration
                if options.rot_inv
                    [cnet, tr_loss(end + 1)] = cnet.rotTrainingIteration(cnet.actvtClass(x_train{1}),cnet.actvtClass(y_train{1}),curr_weights{1});
                else %standard training
                    [cnet, tr_loss(end + 1)] = cnet.trainingIteration(cnet.actvtClass(x_train{1}),cnet.actvtClass(y_train{1}),curr_weights{1});
                end
                
                %training progress printing
                if mod(total_iter, options.display) == 0
                    fprintf('[%s] Iteration %d: lr %.2e, Training loss %.2f\n',datestr(now),total_iter,cnet.optimizer.learningRate,sum(tr_loss)/length(tr_loss));
                    tr_loss = [];
                end
                
                %snapshot saving
                if ~isempty(options.snapshot) && (options.snapshot > 0) && mod(total_iter, options.snapshot) == 0
                    saveCnet(cnet, options, stacksT, stacksV, total_iter, val_loss);
                end
                
                %training abortion due to max iterations reached
                if total_iter == options.max_iter
                    fprintf('[%s] Maximum iteration limit reached.\n',datestr(now));
                    fprintf('[%s] Saving final output.\n',datestr(now));
                    saveCnet(cnet, options, stacksT, stacksV, total_iter, val_loss);
                    cnet = cnet.setParamsToActvtClass(actvtClass);
                    return;
                end
                
                %learning rate adaption
                switch options.lr_policy
                    case 'constant'
                    case 'step'
                        if mod(total_iter,options.step_size) == 0
                            cnet.optimizer.learningRate = cnet.optimizer.learningRate.*options.gamma;
                        end
                    case 'exp'
                        cnet.optimizer.learningRate = base_lr*options.gamma^(total_iter);
                    case 'poly'
                        cnet.optimizer.learningRate = base_lr*(1 - total_iter/options.max_iter)^options.power;
                    otherwise
                        error('lr_policy not implemented');
                end
                total_iter = total_iter + 1;
                skipVal = false;
            end
        end
    end
    if isinf(options.max_iter)
        total_iter = total_iter - 1; %correct for last iteration
        break;
    end
end

%save final cnet
saveCnet(cnet, options, stacksT, stacksV, total_iter, val_loss);

%cast to initial actvt class
cnet = cnet.setParamsToActvtClass(actvtClass);
end

%saving of cnet
function saveCnet(cnet, options, stacksT, stacksV, total_iter, val_loss)
if isfield(options,'snapshot_name') && ~isempty(options.snapshot_name) && options.snapshot > 0
    fprintf('[%s] Iteration %d: Saving snapshot to %s.\n', datestr(now), total_iter, [options.snapshot_name, '_Iter', num2str(total_iter), '.mat']);
    cnet_snapshot = cnet;
    options.data_pre = func2str(options.data_pre);
    options.saveTime = datestr(now);
    cnet_snapshot = cnet_snapshot.setParamsToActvtClass(@single);
    m = matfile([options.snapshot_name, '_Iter', num2str(total_iter), '.mat'],'Writable',true);
    m.cnet = cnet_snapshot;
    m.options = options;
    m.stacksT = stacksT;
    m.stacksV = stacksV;
    if options.val_iter > 0
        m.validationLoss = val_loss;
    end
end
end

%generic loading function
function [raw, target, targetMask] = loadData(path)
m = matfile(path);
raw = single(m.raw);
target = single(m.target);
try
    targetMask = m.targetMask;
catch
    targetMask = false(size(target));
end
end
