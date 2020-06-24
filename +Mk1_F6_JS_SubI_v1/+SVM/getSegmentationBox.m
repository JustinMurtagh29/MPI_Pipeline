function [seg, synCom, vcCom, miCom] = getSegmentationBox(p,segSVM, i,j,k)
% make segmentation to upload in wk
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>
%       : Sahil Loomba <sahil.loomba@brain.mpg.de>

m = load(fullfile(p.local(i,j,k).saveFolder, 'svmPredictions.mat'),...
             'pred17', 'pred17_2', 'pred17_4', 'bbox');
bbox = m.bbox;
m1_pred = m.pred17_2; % use for mitos
m2_pred = m.pred17_4; % use for vesicles

bbox_wk = Util.convertMatlabToWKBbox(bbox); 

predOut = zeros([bbox_wk(4:end), 3]);
predOut(:,:,:,2) = m2_pred(:,:,:,3);
predOut(:,:,:,3) = m1_pred(:,:,:,4);
minArea = 3000;
[seg, synCom, vcCom, miCom] = ...
    Paper.SynEM.MethComp.preprocessSVMForAnnotation(predOut, 0.75, minArea, ...
    [5, 5, 3]);

% further delete small components
for i = 1:3
    tmp = seg == i;
    tmp = regionprops(tmp > 0, 'Area', 'PixelIdxList');
    tmp([tmp.Area] >= minArea) = [];
    idx = cell2mat({tmp.PixelIdxList}');
    seg(idx) = 0;
end

saveSegDataGlobal(segSVM, p.local(i,j,k).bboxSmall(:,1)', seg); 
Util.save([p.local(i,j,k).saveFolder 'svmSegData.mat'], seg,synCom, vcCom, miCom);
end