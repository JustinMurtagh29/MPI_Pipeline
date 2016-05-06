function f = quickPred( cnet, raw, rotPred )
%QUICKPRED Predict raw data on GABA
% INPUT cnet: A Codat.CNN object.
%       raw: (Optonal) 3D array of uint8 for which the prediction is
%           calculated.
%           (Default: File 130.mat of the cortext test set).
%       rotPred: (Optional) Flag to set rotation invariant prediction.
%           (Default: false).
% OUTPUT f: Function handle that will wait for the job to finish and fetch
%           its output. Call pred = f(); to get the job output.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

%prepare input
if ~exist('raw','var') || isempty(raw)
    fprintf('[%s] Loading 130.mat.\n',datestr(now));
    m = load('/u/bstaffle/data/cortexTrainingData/testSet/targetKLEE/130.mat');
    raw = m.raw;
end
raw = (single(raw) - 122)./22;

%job submission
cluster = getCluster('gpu');
cnet = cnet.setConvMode('fft2');
if exist('rotPred','var') && ~isempty(rotPred)
    options.rotPred = rotPred;
end
options.gpuDev = 1;
options.target_size = [50 50 50]; %avoid gpu fft bug
job = createJobPwd(cluster);
createTask(job,@Codat.CNN.Misc.predictCube,1,{cnet, raw, options});
submit(job);

%get output
    function pred = collectOutputs()
        fprintf('[%s] Waiting for job output.\n', datestr(now));
        wait(job);
        fprintf('[%s] Job finished ... fetching outputs.\n',datestr(now));
        out = fetchOutputs(job);
        pred = out{1};
        delete(job);
    end

f = @collectOutputs;
end


