%% Script used for downsampling the datasets
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>


if ispc
    stack = IO.readTiffStack(['E:\workspace\data\backed\Synapse detection\' ...
        'KnottHighResComparison\volumedata.tif']);
    load('E:\workspace\allParametersFileServer.mat')
    stack = stack(:,:,1:500); % local testing
else
    stack = IO.readTiffStack(['/gaba/u/bstaffle/data/Knott_data/' ...
        'volumedata.tif']);
    load('allParameter20141007.mat');
end
raw = Seg.IO.loadRaw(p, 67);
stackD = SynEM.HLRC.downsampleStack(stack, raw, 8./[5, 5, 5]);

if ~ispc
    v = 'v10';
    save('sync/stackD.mat', 'stackD', 'v');
end

%% make knossos hierarchy
IO.writeTif(stackD, 'E:\workspace\data\local\Knott\knott_bad_v10\tifs\', ...
    true);
startFast %script from Flo