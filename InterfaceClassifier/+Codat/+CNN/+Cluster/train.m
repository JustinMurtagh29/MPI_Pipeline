function jobs = train( type, cnets, saveFolder, snapFiles, cluster, comment, varargin )
%TRAIN Cnet training wrapper for cluster.
% INPUT type: Type of cnet to train (see switch below).
%       cnets: Cell array of cnets to train on task.
%       saveFolder: Folder where the snapshot files are save (without 
%                   filesep at the end).
%       snapFiles: Cell array containing the file locations for the cnet
%                  snapshots. If only one name is provided the names will
%                  automatically get index of the cnet attached.
%       cluster: (Optional) A parcluster object.
%                (Default: getCluster('gpu'))
%       comment: (Optional) Variable that will be attached to the options
%                cnn option struct and can contain additional information
%                about the training.
%                (If not specified or empty: No comment will be saved).
%       varargin: Name-value pairs overwriting any of the options listed
%           below. Name should be the fieldname of options.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

%parse input
if ~exist('cluster','var') || isempty(cluster)
    cluster = getCluster('gpu');
end
if exist('comment','var') && ~isempty(comment)
    options.comment = comment;
end

%only important for MMV training
mmvStacks = [1:6 10 13 15 18:20];
mmvStacksWExcl = [1:6 8 10 12 15:17];

switch type
    case 'membrane'
        options.data_pre = func2str(@Codat.CNN.Segmentation.loadData);
        %training data paths
        load('cortexTrainingDataParameter.mat');
        stacks = stacks(~exclude);
        stacksT = {stacks(:).stackFile};
        %test data path
        load('cortexTestDataParameter.mat')
        stacks = stacks(~exclude);
        stacksV = {stacks(1:10).targetFile};
        
    case 'membraneWeighted'
        options.data_pre = func2str(@(x)Codat.CNN.Segmentation.loadDataWeighted(x,8));
        %training data paths
        load('cortexTrainingDataParameter.mat');
        stacks = stacks(~exclude);
        stacksT = {stacks(:).stackFile};
        %test data path
        load('cortexTestDataParameter.mat')
        stacks = stacks(~exclude);
        stacksV = {stacks(1:10).targetFile};
        
	case 'membraneWeightedsT'
        options.data_pre = func2str(@(x)Codat.CNN.Segmentation.loadDataWeighted(x,10,true));
        %training data paths
        load('cortexTrainingDataParameter.mat');
        stacks = stacks(~exclude);
        stacksT = {stacks(:).stackFile};
        %test data path
        load('cortexTestDataParameter.mat')
        stacks = stacks(~exclude);
        stacksV = {stacks(1:10).targetFile};
        
    case 'membraneRMMVWeighed'
        options.data_pre = func2str(@(x)Codat.CNN.Segmentation.loadDataRMMV(x,8));
        %training data paths
        load('cortexTrainingDataParameterMMV.mat');
        stacks = stacks(~exclude);
        stacksT = {stacks(:).stackFile};
        %test data path
        load('cortexTestDataParameterMMV.mat')
        stacks = stacks(~exclude);
        stacksV = {stacks(1:10).targetFile};
        
    case 'membraneRMMVWeighedsT'
        options.data_pre = func2str(@(x)Codat.CNN.Segmentation.loadDataRMMV(x,10,true));
        %training data paths
        load('cortexTrainingDataParameterMMV.mat');
        stacks = stacks(~exclude);
        stacksT = {stacks(:).stackFile};
        %test data path
        load('cortexTestDataParameterMMV.mat')
        stacks = stacks(~exclude);
        stacksV = {stacks(1:10).targetFile};
        
    case 'membraneSig'
        options.data_pre = func2str(@Codat.CNN.Segmentation.loadDataSig);
        %training data paths
        load('cortexTrainingDataParameter.mat');
        stacks = stacks(~exclude);
        stacksT = {stacks(:).stackFile};
        %test data path
        load('cortexTestDataParameter.mat')
        stacks = stacks(~exclude);
        stacksV = {stacks(1:10).targetFile};
        
    case 'membraneWeightedSig'
        options.data_pre = func2str(@(x)Codat.CNN.Segmentation.loadDataSigWeighted(x,8));
        %training data paths
        load('cortexTrainingDataParameter.mat');
        stacks = stacks(~exclude);
        stacksT = {stacks(:).stackFile};
        %test data path
        load('cortexTestDataParameter.mat')
        stacks = stacks(~exclude);
        stacksV = {stacks(1:10).targetFile};

    case 'MMV'
        options.data_pre = func2str(@Codat.CNN.Segmentation.loadDataMT);
        %training & test paths
        load('cortexTrainingDataParameter.mat');
        stacksT = {stacks(mmvStacks).stackFile};
        stacksV = {stacks(10).stackFile};
        
    case 'MMVWeighted'
        options.data_pre = func2str(@(x)Codat.CNN.Segmentation.loadDataMMVWeighted(x,8));
        %training & test paths
        load('cortexTrainingDataParameter.mat');
        stacksT = {stacks(mmvStacks).stackFile};
        stacksV = {stacks(10).stackFile};
        
    case 'MMVMasked'
        options.data_pre = func2str(@Codat.CNN.Segmentation.loadDataMTMasked);
        %training data paths
        load('cortexTrainingDataParameter.mat');
        stacks = stacks(~exclude);
        stacksT = [stacks(:).stackFile repmat({stacks(mmvStacksWExcl).stackFile},1,ceil(length(stacks)/length(mmvStacksWExcl)) - 1)];
        stacksV = {stacks(mmvStacksWExcl).stackFile};
        %test data paths
        load('cortexTestDataParameter.mat')
        stacks = stacks(~exclude);
        stacksV = [stacksV {stacks(1:10).targetFile}];

    case 'MMVMaskWeigh'
        options.data_pre = func2str(@(x)Codat.CNN.Segmentation.loadDataMMVMaskWeigh(x,8));
        %training data paths
        load('cortexTrainingDataParameter.mat');
        stacks = stacks(~exclude);
        stacksT = [stacks(:).stackFile repmat({stacks(mmvStacksWExcl).stackFile},1,ceil(length(stacks)/length(mmvStacksWExcl)) - 1)];
        stacksV = {stacks(mmvStacksWExcl).stackFile};
        %test data paths
        load('cortexTestDataParameter.mat')
        stacks = stacks(~exclude);
        stacksV = [stacksV {stacks(1:10).targetFile}];
        
    case 'MMVMaskWeighSig'
        options.data_pre = func2str(@(x)Codat.CNN.Segmentation.loadDataMMVMaskWeighSig(x,8));
        %training data paths
        load('cortexTrainingDataParameter.mat');
        stacks = stacks(~exclude);
        stacksT = [stacks(:).stackFile repmat({stacks(mmvStacksWExcl).stackFile},1,ceil(length(stacks)/length(mmvStacksWExcl)) - 1)];
        stacksV = {stacks(mmvStacksWExcl).stackFile};
        %test data paths
        load('cortexTestDataParameter.mat')
        stacks = stacks(~exclude);
        stacksV = [stacksV {stacks(1:10).targetFile}];
        
    case 'MMVWeightedSig'
        options.data_pre = func2str(@(x)Codat.CNN.Segmentation.loadDataMMVWeighSig(x,8));
        %training & test paths
        load('cortexTrainingDataParameter.mat');
        stacksT = {stacks(mmvStacks).stackFile};
        stacksV = {stacks(10).stackFile};
        
    case 'Synapses'
        options.data_pre = func2str(@Codat.CNN.Synapses.loadDataWeighted);
        load('synapseVoxelTrainingDataParameter.mat');
        stacksT = {stacks(:).stackFile};
        stacksV = {stacks([12 16 17 21 35 36 44 73 79 96 141 157]).stackFile};
        
    otherwise
        error('Unknown type %s.\n',type);
end

%define training options
options.display = 125;
options.tr_size = [20, 20, 20];
options.val_size = [50, 50, 50];
options.val_iter = 2500;
options.val_fwd_alg = 'fft2';
options.gpuDev = true;
options.max_iter = 200000;
options.augment = false;
options.snapshot = 10000;
options.save_imp = false;
options.class_ratio = [];
options.skip_single_class_cubes = false;
if isa(cnets{1}.optimizer,'Codat.Optimizer.adam')
    options.lr_policy = 'constant';
else
    options.lr_policy = 'step';
end
options.step_size = 50000;
options.gamma = 0.1;
options.rot_inv = false;

%parse name/value inputs to change standard options
nArgs = length(varargin);
optionNames = fieldnames(options);
if floor(nArgs/2) ~= nArgs/2
    error('Optional arguments have to be name/value pairs');
end
for pair = reshape(varargin,2,[])
    name = pair{1};
    if any(strcmp(name,optionNames))
        options.(name) = pair{2};
    else
        error('%s is not a recognized parameter name',name);
    end
end

%attach number of cnet in case of a single name
if ischar(snapFiles)
    tmp = snapFiles;
    snapFiles = cell(length(cnets),1);
    for i = 1:length(cnets)
        snapFiles{i} = [tmp '_' num2str(i)];
    end
elseif length(snapFiles) == 1 && length(cnets) > 1
    for i = 1:length(cnets)
        snapFiles{i} = [snapFiles{1} '_' num2str(i)];
    end
end

%submit jobs
for i = 1:length(cnets)
    cnet = cnets{i};
    options.snapshot_name = [saveFolder filesep snapFiles{i}];
    jobs(i) = createJobPwd(cluster);
    jobs(i).Name = [snapFiles{i} '_' num2str(jobs(i).ID)];
    t = createTask(jobs(i),@cnet.train,2,{stacksT,stacksV,options});
    t.CaptureDiary = true;
    submit(jobs(i));
end


end

