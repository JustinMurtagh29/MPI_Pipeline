function f = validateCnets( cnets, type, cluster )
%VALIDATECNETS Validate multiple cnets.
% INPUT cnets: [Nx1] cell array where each cell contains a Codat.CNN.cnn
%           object or path to folder where cnets are stored (in this case
%           all m-files in the folder are loaded and each should contain a
%           Codat.CNN.cnn object saved with variable name "cnet").
%       type: Type of cnet (see switch below).
%       cluster: (Optional) Parcluster object.
%                (Default: getCluster('gpu'))
% NOTE Results are saved in folder with filename loss.mat.
%
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if ischar(cnets)
    folder = cnets;
    s = what(folder);
    cnets = cell(length(s.mat),1);
    for i = 1:length(s.mat)
        m = matfile([folder filesep s.mat{i}]);
        cnets{i} = m.cnet;
    end
end

options.gpuDev = true;
options.val_fwd_alg = 'fft2';
options.val_size = [50 50 50];
switch type
    case 'membrane'
        m = load('/gaba/u/bstaffle/code/workspace/cortexTestDataParameter.mat');
        stacks = m.stacks(~m.exclude);
        stacksV = {stacks(1:10).targetFile};
        options.data_pre = func2str(@(x)Codat.CNN.Segmentation.loadDataWeighted(x,8));
    case 'membraneSig'
        m = load('/gaba/u/bstaffle/code/workspace/cortexTestDataParameter.mat');
        stacks = m.stacks(~m.exclude);
        stacksV = {stacks(1:10).targetFile};
        options.data_pre = func2str(@(x)Codat.CNN.Segmentation.loadDataSigWeighted(x,8));
    case 'MMV'
        m = load('/gaba/u/bstaffle/code/workspace/cortexTrainingDataParameter.mat');
        stacksV = {m.stacks([1:6 10]).stackFile};
        options.data_pre = func2str(@(x)Codat.CNN.Segmentation.loadDataMMVWeighted(x,8));
    case 'MMVSig'
        m = load('/gaba/u/bstaffle/code/workspace/cortexTrainingDataParameter.mat');
        stacksV = {m.stacks([1:6 10]).stackFile};
        options.data_pre = func2str(@(x)Codat.CNN.Segmentation.loadDataMMVWeighSig(x,8));
    case 'MMVMasked'
        m = load('/gaba/u/bstaffle/code/workspace/cortexTestDataParameter.mat');
        stacks = m.stacks(~m.exclude);
        stacksV = {stacks(1:10).targetFile};
        options.data_pre = func2str(@(x)Codat.CNN.Segmentation.loadDataMMVMaskWeigh(x,8));
    case 'MMVMaskedSig'
        m = load('/gaba/u/bstaffle/code/workspace/cortexTestDataParameter.mat');
        stacks = m.stacks(~m.exclude);
        stacksV = {stacks(1:10).targetFile};
        options.data_pre = func2str(@(x)Codat.CNN.Segmentation.loadDataMMVMaskWeighSig(x,8));
    case 'Synapses'
        m = load('synapseVoxelTrainingDataParameter.mat');
        stacksV = {m.stacks([12 16 17 21 35 36 44 73 79 96 141 157]).stackFile};
        options.data_pre = func2str(@Codat.CNN.Synapses.loadDataMasked);
    otherwise
        error('Unknown type %s.\n',type);
end

if ~exist('cluster','var') || isempty(cluster)
    cluster = getCluster('gpu');
end

fprintf('[%s] Submitting jobs to cluster.\n', datestr(now));
inputCell = cell(length(cnets),1);
for i = 1:length(inputCell)
    inputCell{i} = {cnets{i}, stacksV, options};
end
job = startJob(cluster,@validateWrapper, inputCell, 1);

    function t = collectOutputs()
        fprintf('[%s] Waiting for job output.\n', datestr(now));
        wait(job);
        try
            out = fetchOutputs(job);
            loss = cell2mat(out);
        catch
            warning('Some error occured during loading. Loading tasks separately');
            tasks = job.Tasks;
            loss = zeros(length(tasks),1);
            for j = 1:length(tasks)
                try
                    loss(j) = tasks(j).OutputArguments{1};
                catch
                    warning('Error occured during task %d.', j);
                end
            end
        end
        
        t = loss;
        if exist('folder','var')
            files = s.mat;
            t = table(loss,files);
            t = sortrows(t, 'loss', 'ascend');
            targetFile = [folder filesep 'valLoss.mat'];
            fprintf('[%s] Saving results to %s.\n', datestr(now), targetFile);
            m = matfile(targetFile, 'Writable', true);
            m.stacksV = stacksV;
            m.loss = t;
        end
    end

f = @collectOutputs;

end

function loss = validateWrapper(cnet, stacksV, options)
if options.gpuDev
    cnet = cnet.setParamsToActvtClass(@gpuArray);
end
loss = cnet.validate(stacksV,options);
end
