function result = speedTest( cnet, cluster, targetSize, mode, profileDir )
%SPEEDTEST Test speed of forward pass and backprop.
% INPUT cnet: A Codat.CNN object.
%       cluster: A matlab parcluster object.
%       targetSize: [3x1] numerical array specifying the target size in
%           each dimension.
%       mode: String specifying the calculation device. Choices are 'cpu',
%           'gpu' or 'mixed' for forward pass on GPU and backward pass on
%           CPU.
%       profileDir: Directory to save profile. If not specified or empty no
%                   profile will be made.
% OUTPUT result: Struct with speed test results.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if ~exist('profileDir','var')
    profileDir = [];
end

input = randn([targetSize(1:3) + cnet.border,cnet.featureMaps(1)],'single');
job = createJobPwd(cluster);
switch mode
    case 'cpu'
        createTask(job,@CPUSpeedTest,1,{cnet, input, profileDir});
    case 'gpu'
        createTask(job,@GPUSpeedTest,1,{cnet, input, profileDir});
    case 'mixed'
        createTask(job,@mixedSpeedTest,1,{cnet, input, profileDir});
end
submit(job);
wait(job);
out = fetchOutputs(job);
result = out{1};

end

function result = CPUSpeedTest(cnet, input, profileDir)

cnet.isTraining = true;

if any(strcmp(cnet.convAlg,{'fft1','fft2'})) %fft parameter search
    [ activity, dropoutMask, mpInd ] = forwardPass( cnet, input );
    target = cnet.actvtClass(randn(size(activity{end})));
    targetMask = false(size(target));
    backprop( cnet, activity, dropoutMask, mpInd, target,targetMask );
end

if ~isempty(profileDir)
    profile on
end

tic
[ activity, dropoutMask, mpInd ] = forwardPass( cnet, input );
result.t_fwd = toc;

target = cnet.actvtClass(randn(size(activity{end})));
targetMask = false(size(target));

tic
backprop( cnet, activity, dropoutMask, mpInd, target, targetMask );
result.t_bwd = toc;

if ~isempty(profileDir)
    p = profile('info');
    profsave(p,profileDir);
end
end

function result = GPUSpeedTest(cnet, input, profileDir)

cnet.isTraining = true;
input = gpuArray(input);
cnet = cnet.setParamsToActvtClass(@gpuArray);

if any(strcmp(cnet.convAlg,{'fft1','fft2'})) %fft parameter search
    [ activity, dropoutMask, mpInd ] = forwardPass( cnet, input );
    target = cnet.actvtClass(randn(size(activity{end})));
    targetMask = false(size(target));
    backprop( cnet, activity, dropoutMask, mpInd, target, targetMask );
end

if ~isempty(profileDir)
    profile on
end

tic
[ activity, dropoutMask, mpInd ] = forwardPass( cnet, input );
result.t_fwd = toc;

target = gpuArray.randn(size(activity{end}),'single');
targetMask = gpuArray.false(size(target));

tic
backprop( cnet, activity, dropoutMask, mpInd, target, targetMask );
result.t_bwd = toc;

if ~isempty(profileDir)
    p = profile('info');
    profsave(p,profileDir);
end
end

function result = mixedSpeedTest(cnet, input, profileDir)

cnet.isTraining = true;
input = gpuArray(input);
cnet = cnet.setParamsToActvtClass(@gpuArray);

if any(strcmp(cnet.convAlg,{'fft1','fft2'})) %fft parameter search
    forwardPass( cnet, input );
end

if ~isempty(profileDir)
    profile on
end

tic
[ activity, dropoutMask, mpInd ] = forwardPass( cnet, input );
result.t_fwd = toc;

target = randn(size(activity{end}),'single');
targetMask = false(size(target));

tic
activity = cellfun(@(x)gather(x),activity,'UniformOutput',false);
mpInd = cellfun(@(x)gather(x),mpInd,'UniformOutput',false);
dropoutMask = cellfun(@(x)gather(x),dropoutMask,'UniformOutput',false);
cnet = cnet.setParamsToActvtClass(@single);
backprop( cnet, activity, dropoutMask, mpInd, target, targetMask );
result.t_bwd = toc;

if ~isempty(profileDir)
    p = profile('info');
    profsave(p,profileDir);
end
end

