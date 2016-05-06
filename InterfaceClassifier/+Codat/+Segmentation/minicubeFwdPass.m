function job = minicubeFwdPass( parameter, cluster )
% classification of subergions within the data set
% (for segmentation optimization and GP training)
% INPUT parameter: A segmentation parameter struct (see also
%           setParameters07x2).
%       cluster: (Optional) A matlab cluster object.
%           (Default: getCluster('cpu') or getCluster('gpu') depending on
%               parameter.cnn.GPU)
% Author: Manuel Berning <manuel.berning@brain.mpg.de>
% Modified by Benedikt Staffler <benedikt.staffler@brain.mpg.de>
    
% Settings for Benedikt's CNN
options.gpuDev = parameter.cnn.GPU;
options.convAlg = 'fft2';
options.target_size = [128 128 128];

% Initalize counter variable (each KNOSSOS cube will be one task)
nrTasks = 0;
for tr=1:length(parameter.local)
    bbox = parameter.local(tr).bboxBig;
    % Bbox inidcating lower,left,front corner of Knossos cubes
    bbox = floor((bbox-1)./128).*128+1;
    for x=bbox(1,1):128:bbox(1,2)
        for y=bbox(2,1):128:bbox(2,2)
            for z=bbox(3,1):128:bbox(3,2)
                % Bbox of this Knossos cube
                thisBbox = [x x+127; y y+127; z z+127];
                nrTasks = nrTasks + 1;
                inputCell{nrTasks} = {parameter.cnn.saveFile, thisBbox, options, parameter.raw, parameter.local(tr).class};
            end
        end
    end
end
functionH = @jobWrapper;

if ~exist('cluster','var') || isempty(cluster)
    if parameter.cnn.GPU
        cluster = getCluster('gpu');
    else
        cluster = getCluster('cpu');
    end
end
job = startJob(cluster, functionH, inputCell,[],'classification');

end

function jobWrapper(cnnFile, ROI, options, knossosConfRaw, knossoConfClass)
m = load(cnnFile);
cnet = m.cnet;
pred = Codat.CNN.Misc.predictROI(cnet, ROI, options, knossosConfRaw);
writeKnossosRoi(knossoConfClass.root, knossoConfClass.prefix, ROI(:,1)', pred, 'single','','noRead');
end

