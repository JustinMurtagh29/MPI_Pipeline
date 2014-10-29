clc;
cd /home/mberning/code/active/;

%% load training data
% weights = cell(9,1);
% edges = cell(9,1);
% labelEdges = cell(9,1);
% load dataWithGroundTruth.mat;
% raw = single(raw);
% raw = raw - mean(raw(:));
% raw = raw ./ std(raw(:));
% for x=1:3
%     for y=1:3
%         display([num2str(1+(x-1)*100) ' ' num2str(x*100) ' ' num2str(1+(y-1)*100) ' ' num2str(y*100)]);
%         tic;
%         affT = aff(1+(x-1)*100:x*100,1+(y-1)*100:y*100,1:100);
%         rawT = raw(1+(x-1)*100:x*100,1+(y-1)*100:y*100,1:100);
%         seg_manuellT = seg_manuell(1+(x-1)*100:x*100,1+(y-1)*100:y*100,1:100);
%         clear seg;
%         %%calculate edges, label for edges and weights
%         seg = redoSegmentation( affT );
%         trueEdges = approximateGroundTruth(seg, seg_manuellT);
%         edges{sub2ind([3 3],x,y)} = findEdges(seg);
%         labelEdges{sub2ind([3 3],x,y)} = findLabels( edges{sub2ind([3 3],x,y)}, trueEdges );
%         [weights{sub2ind([3 3],x,y)}, ~] = featureDesign(rawT, affT, seg, edges{sub2ind([3 3],x,y)});
%         time = toc;
%         display([num2str(ceil(time/60)) ' min']);
%     end
% end

%%
load forMoritz.mat;
aff = aff(1:100,1:100,1:100);
raw = single(raw(1:100,1:100,1:100));
raw = raw - mean(raw(:));
raw = raw ./ std(raw(:));

%%
profile on;
display('large test cube');
clear seg; 
seg = redoSegmentation( aff );
edges2 = findEdges(seg);
[weights2, weightLabels] = featureDesign(raw, aff, seg, edges2);
save('wipSpeedTest.mat');
profile viewer
profile off;

%%
profile on;
result = calcFeatures(aff);
profile viewer
profile off;