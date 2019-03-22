function [seg, synCom, vcCom, miCom] = getSegmentationBox(p, cubeId)
% make segmentation to upload in wk
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>
%       : Sahil Loomba <sahil.loomba@brain.mpg.de>

m = load(fullfile(p.local(cubeId).saveFolder, 'svmPredictions.mat'),...
             pred17, pred17_2, pred17_4, bbox);
m1 = m.pred17_2; % use for mitos
m2 = m.pred17_4; % use for vesicles

bbox_wk = Util.convertMatlabToWKBbox(bbox); 

pred = zeros([bbox_wk(4:end), 3]);
pred(:,:,:,2) = m2.pred{i}(:,:,:,3);
pred(:,:,:,3) = m1.pred{i}(:,:,:,4);
minArea = 3000;
[seg, synCom, vcCom, miCom] = ...
    Paper.SynEM.MethComp.preprocessSVMForAnnotation(pred, 0.75, minArea, ...
    [5, 5, 3]);

% further delete small components
for i = 1:3
    tmp = seg == i;
    tmp = regionprops(tmp > 0, 'Area', 'PixelIdxList');
    tmp([tmp.Area] >= minArea) = [];
    idx = cell2mat({tmp.PixelIdxList}');
    seg(idx) = 0;
end

seg_bbox = seg;

end
