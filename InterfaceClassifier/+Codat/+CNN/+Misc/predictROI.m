function [pred, raw] = predictROI( cnet, ROI, options, knossosConf )
%PREDICTROI Prediction for ROI.
% INPUT cnet: A Codat.CNN object.
%       ROI: Region of interest as 3x2 coordinates for prediction
%            (e.g. [500 1000; 500 1000;500 750]).
%            Prediction will be made for the whole ROI, i.e. an additional
%            boundary of cnet.border/2 in each dimension of the ROI is
%            loaded
%       options: Options struct with field
%       	gpuDev: see cnet.train
%       	target_size: 3x1 array.
%               Input cube is tiled such that each forward
%               pass produces an output of target_size.
%               Targets for one cube are stiched together
%               afterwards. (Default is the whole ROI).
%           convAlg (Default is current convAlg of cnet)
%           f_norm: (Optional) Function handle for normalization of
%               raw data.
%               (Default: (raw - 122)./22)
%       knossosConf: Struct containing the fields
%           'root': Path to knosso hierarchy root folder
%           'prefix': Cube fileprefix
%           (e.g. use p.raw for a segmentatin parameter struct p)
% OUTPUT pred: Prediction for ROI.
%        raw: ROI raw data (uint8) including boundary.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

%set options
if options.gpuDev
    % Comment MB: not needed anymore after GPU exclusive mode
    % chooseGPUDev();
    cnet = cnet.setParamsToActvtClass(@gpuArray);
end
if ~isfield(options,'target_size') || isempty(options.target_size)
    options.target_size = (ROI(:,2) - ROI(:,1) + 1)';
elseif isrow(options.target_size)
    options.target_size = options.target_size';
end
if isfield(options,'convAlg')
    cnet = cnet.setConvMode(options.convAlg);
end
if ~isfield(options,'f_norm')
    options.f_norm = @(x)(single(x) - 122)./22;
end


%tile ROI
[y,x,z] = meshgrid(ROI(2,1):options.target_size(1):ROI(2,2), ...
                   ROI(1,1):options.target_size(2):ROI(1,2), ...
                   ROI(3,1):options.target_size(3):ROI(3,2));
coords = [x(:),y(:),z(:)]';
pred = cell(size(x));

%predict each tile
for tile = 1:numel(x)
    %add cnet border to tile ROI
    tileROI = [coords(:,tile) - cnet.border'./2, ...
        min(coords(:,tile) + options.target_size + cnet.border'./2 - 1,ROI(:,2) + cnet.border'./2)];
    raw = readKnossosRoi(knossosConf.root,knossosConf.prefix,tileROI,'uint8');
    raw = options.f_norm(raw);
    if options.gpuDev
        raw = gpuArray(raw);
    end
    pred{tile}= gather(cnet.predict(raw));
end
pred = cell2mat(pred);

%load full input cube if required
if nargout == 2
    raw = readKnossosRoi(knossosConf.root,knossosConf.prefix,ROI + [-cnet.border'/2, cnet.border'/2],'uint8');
end

end
