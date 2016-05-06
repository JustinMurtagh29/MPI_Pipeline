function job = extendCortexTrainingData(stacks, cnet, options, cluster)
%EXTENDCORTEXTRAININGDATA Extend the cortex training data by the output of the
% provided cnet.
% INPUT stacks: Cortex training data parameter struct.
%       cnet: A Codat.CNN object.
%       options: Structure which must/can have the following fields
%           'outputName': (Optional) String specifying the name of the
%           	additional variable saved into each cortex training data
%           	mat-file.
%               (Default: 'feat')
%           'outputSize': (Optional) [3x1] numerical array specifying the
%           	size to which the output of the cnet prediction is cropped
%               before saving it.
%               (Default or []: Prediction is not cropped).
%           'gpuDev': (Optional) Bool specifying whether to use gpu.
%           	(Default: true).
%           'convAlg': (Optional) see Codat.CNN.cnn convAlg
%               (Default: 'fft2')
%       cluster: (Optional) A matlab parcluster object.
%                (Default or []: getCluster('gpu') or getCluster('cpu')
%                 depending on options.gpuDev)
% OUTPUT job: The resulting job object.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

%parse inputs
if ~exist('options','var') || isempty(options)
    options = struct;
end
if ~isfield(options,'outputName') || isempty(options.outputName)
    options.outputName = 'feat';
end
if isfield(options,'convAlg') && ~isempty(options.convAlg)
    cnet = cnet.setConvMode(options.convAlg);
else
    cnet = cnet.setConvMode('fft2');
end
if ~isfield(options,'outputSize') || isempty(options.outputSize)
    options.outputSize = [];
else
    predSize = [200, 200, 150] - cnet.border;
    border = predSize - options.outputSize;
    if any(border < 0)
        error('cnet border is too large to achieve desired output size');
    elseif any((border./2) ~= floor(border./2))
        error('cnet output does not allow for symmetric cropping of prediction.')
    end
end
if ~isfield(options,'gpuDev') || isempty(options.gpuDev)
    options.gpuDev = true;
end

%get cluster
if (~exist('cluster','var') || isempty(cluster)) && options.gpuDev
    cluster = getCluster('gpu');
elseif (~exist('cluster','var') || isempty(cluster)) && ~options.gpuDev
    cluster = getCluster('cpu');
end

%create and submit job
inputCell = cell(length(stacks),1);
for i = 1:length(stacks)
    inputCell{i} = {stacks(i).stackFile, cnet, options};
end
job = startJob(cluster,@jobWrapper,inputCell,0);
end

function jobWrapper(file, cnet, options)
%Wrapper for job submitted to cluster. Uses cnet to predict from
% cortextTrainingData stack raw file and saves the output.

%calculate cnet prediction
m = matfile(file,'Writable',true);
raw = (single(m.raw) - 122)./22;
if options.gpuDev
    cnet = cnet.setParamsToActvtClass(@gpuArray);
    raw = gpuArray(raw);
end
pred = gather(cnet.predict(raw));

%crop prediction if necessary
if ~isempty(options.outputSize)
    b = (size(pred) - options.outputSize)./2;
    pred = pred(b(1) + 1:end - b(1),b(2) + 1:end - b(2),b(3) + 1:end - b(3),:);
end
%save to cortexTraining stack
m.(options.outputName) = pred;
end
