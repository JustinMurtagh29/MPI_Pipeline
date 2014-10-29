clc;
cd C:\Users\mberning\Desktop\active;
load dataWithGroundTruth.mat;
aff = aff(1:100,1:100,1:100);
raw = raw(1:100,1:100,1:100);
seg = seg(1:100,1:100,1:100);
seg_manuell = seg_manuell(1:100,1:100,1:100);

%% Redo segmentation (to really oversegment)
map = imimposemin(imcomplement(aff), imextendedmin(imcomplement(aff), .18));
seg2 = watershed(map, 26);

%% Find correspondence of manuell and automated segmentation
a = unique(seg_manuell);
a(a == 0) = [];
for i=1:length(a)
    withinSegment{i} = seg2(seg_manuell == a(i));
    withinUnique{i} = unique(withinSegment{i});
    withinUnique{i}(withinUnique{i} == 0) = [];  
    counts{i} = histc(withinSegment{i},withinUnique{i}); 
    withinUnique{i}(counts{i} < 100) = [];
    counts{i}(counts{i} < 100) = [];
end

nrIDs = length(unique(seg2));
nrIDs(nrIDs == 0) = [];
trueEdges = [];
for i=1:size(withinUnique,2)
    for j=1:size(withinUnique{i},1)
        for k=1:size(withinUnique{i},1)
            if k~=j
                trueEdges(end+1,:) = [withinUnique{i}(k) withinUnique{i}(j)];
            end
        end
    end
end
b = cell2mat(withinUnique');
uniqueB = unique(b);
countsB = histc(b,uniqueB);
