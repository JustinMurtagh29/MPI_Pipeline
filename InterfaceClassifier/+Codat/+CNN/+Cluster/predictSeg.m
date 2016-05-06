function job = predictSeg( pFile, cnetFile )
%PREDICTSEG CNN prediction for local segmentation cubes.
% INPUT p: Path to segmentation parameter struct.
%       cnetFile: Path to .mat file containing a Codat.CNN.cnn object saved
%           with name cnet.
% OUTPUT job: Matlab job object.
%
% NOTE This function save the variable cnetFile to p.saveFolder for
%      documenation.
%
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

m = load(pFile);
p = m.p;
save([p.saveFolder, 'MMV_parameter.mat'], 'cnetFile');
inputCell = cell(numel(p.local),1);
for i = 1:numel(p.local);
    inputCell{i} = {pFile, cnetFile, i};
end

cluster = Cluster.getCluster('-pe openmp 1 -l h_vmem=12G,h_rt=100:00:00 -p -500');
job = startJob(cluster,@jobWrapper,inputCell);
end

function jobWrapper(pFile, cnetFile, i)

%load segmentation parameters
m = load(pFile);
p = m.p;

%load cnet
m = load(cnetFile);
cnet = m.cnet;

%load raw
bboxRaw = p.local(i).bboxSmall + [-cnet.border'./2,cnet.border'./2];
raw = readKnossosRoi(p.raw.root,p.raw.prefix,bboxRaw,'uint8');
if all(~raw)
    error('The raw file read from %s is empty.',p.raw.root);
end
raw = (single(raw) - 122)./22;

%prediction
options.gpuDev = false;
options.val_fwd_alg = 'fft2';
options.target_size = [128, 128, 128];
pred = Codat.CNN.Misc.predictCube(cnet, raw, options);

%save pred
m = matfile([p.local(i).saveFolder 'Pred_MMV.mat'],'Writable',true);
m.pred = pred;
end

