function job = predictROI(cnet, ROI, knossosConf, outputFile, options, cluster)
%PREDICTROI Predict ROI on cluster.
% INPUT cnet: A Codat.CNN object.
%       ROI: Region of interest as 3x2 coordinates for prediction
%            (e.g. [500 1000; 500 1000;500 750]).
%            Prediction will be made for the whole ROI, i.e. an additional
%            boundary of cnet.border/2 in each dimension of the ROI is
%            loaded
%       knossosConf: Struct containing the fields
%           'root': Path to knosso hierarchy root folder
%           'prefix': Cube fileprefix
%           (e.g. use p.raw for a segmentatin parameter struct p)
%       outputFile: Full path to output file.
%       options: Options struct with field
%                gpuDev: see cnet.train (Default: true)
%                target_size: 3x1 array.
%                             Input cube is tiled such that each forward
%                             pass produces an output of target_size.
%                             Targets for one cube are stiched together
%                             afterwards. (Default is the whole ROI).
%                val_fwd_alg (Default is current fwdAlg of cnet)
%       cluster: (Optional) A matlab parallel.cluster object.
%           (Default: getCluster('gpu'))
% OUTPUT job: A matlab job object.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if ~exist('options','var') || isempty(options)
    options = struct;
end
if ~isfield(options,'gpuDev') || isempty(options.gpuDev)
    options.gpuDev = true;
end
if ~exist('cluster','var') || isempty(cluster)
    cluster = getCluster('gpu');
end
inputCell = {cnet, ROI, options, knossosConf, outputFile};
job = startJob(cluster, @jobWrapper,inputCell,0);

end

function jobWrapper(cnet, ROI, options, knossosConf, outputFile)
    pred = Codat.CNN.Misc.predictROI(cnet, ROI, options, knossosConf);
    m = matfile(outputFile,'Writable',true);
    m.pred = pred;
end
