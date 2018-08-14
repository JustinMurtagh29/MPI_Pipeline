function detectNucleiinBoxes()
%%
% Detection of nuclei in mag1. The function uses the knowledge of the
% coordinates for a bounding box around each nucleus from the previously
% proceeded detection in mag4:
load('/gaba/u/scchr/MouseP28/NucleiVesselDetection/NucleiCoordinates.mat','rp')

addpath(genpath('/u/scchr/Repositories/auxiliaryMethods/'));
addpath(genpath('/u/scchr/Repositories/pipeline/'));
addpath(genpath('/u/scchr/Repositories/nuclear_pores/'));

prefix = '2012-11-23_ex144_st08x2New_mag1';
root = '/u/scchr/wKcubes/2012-11-23_ex144_st08x2New/color/1/';

% Load mask and vessel data. The vessel need to be excluded during
% processing.
load('/gaba/u/scchr/MouseP28/NucleiVesselDetection/mask.mat');
load('/gaba/u/scchr/MouseP28/NucleiVesselDetection/Vessels_mag4.mat');
Mask = permute(Mask, [1 3 2]);
vessels = permute(vessels, [1 3 2]);
%%
margin = [15, 15, 15]; %margin around bounding box of nucleus

mag1bbox = [170, 440, 32, 8405, 5550, 7745];
mag1bbox = Util.convertWebknossosToMatlabBbox(mag1bbox);

%%
for i=1:length(rp)
    % get bounding box around nucleus with margin for cropped mag4
    bboxM4Cropped = [round(rp(i).BoundingBox(1:3)) - margin; ...
        round(rp(i).BoundingBox(1:3) + rp(i).BoundingBox(4:6)) ...
        + margin]';
    if bboxM4Cropped(1)<1
        bboxM4Cropped(1) = 1;
    end
    if bboxM4Cropped(2)<1
        bboxM4Cropped(2) = 1;
    end
    if bboxM4Cropped(3)<1
        bboxM4Cropped(3) = 1;
    end
    if bboxM4Cropped(4)>2101
        bboxM4Cropped(4) = 2101;
    end
    if bboxM4Cropped(5)>1388
        bboxM4Cropped(5) = 1388;
    end
    if bboxM4Cropped(6)>1937
        bboxM4Cropped(6) = 1937;        
    end
    % corresponding box in mag1, consider bbox of whole dataset (mag1bbox)
    bbox = bsxfun(@plus,bsxfun(@plus,4.*bboxM4Cropped,[-1, 2]), ...
             mag1bbox(:,1));
%     bbox = bsxfun(@plus,4.*bboxM4Cropped, mag1bbox(:,1));
    mask = Util.getBbox(Mask, bboxM4Cropped);
    vessel = Util.getBbox(vessels, bboxM4Cropped);
    BoxList{i,1} = {bbox, mask, vessel, root, prefix, i};
end

%%
% Several bboxes have to be manually modified since they are consisting more
% than one nucleus
empty = [15, 16, 31, 42, 81, 87, 88, 97, 123, 135, 144, 148, 159, 162, 171, 180, 186, 193];

for i=1:length(empty)
    BoxList{empty(i)} = {};
end

box{1} = [270, 1149; 1856, 2487; 128, 591]; % 16
box{2} = [554, 1569; 2378, 3143; 604, 1224]; % 31
box{3} = [1155, 1996; 2459, 3310; 1241, 1613]; % 42
box{4} = [1810, 2590; 2317, 3465; 1180, 1576]; % 42
box{5} = [1733, 2624; 1575, 2456; 1215, 1508]; % 42
box{6} = [3386, 4484; 4289, 5077; 2766, 3234]; % 81
box{7} = [3416, 4252; 4639, 5618; 2620, 3078]; % 81
box{8} = [6586, 7601; 3660, 4363; 3380, 4011]; % 97
box{9} = [1949, 3160; 1730, 2719; 4594, 5070]; % 123
box{10} = [2053, 2958; 2443, 3300; 4785, 5160]; % 123
box{11} = [2963, 3979; 3801, 4678; 5216, 5583]; % 135
box{12} = [2970, 3770; 4498, 5076; 5066, 5351]; % 135
box{13} = [816, 1658; 1264, 2051; 5596,5868]; % 148
box{14} = [5322, 6169; 3408, 3875; 6036, 6488]; % 159
box{15} = [445, 1287; 4679, 5693; 6286, 6553]; % 162
box{16} = [421, 1523; 4384, 5383; 6492, 7003]; % 162
box{17} = [652, 1419; 1362, 2073; 7092, 7346]; % 186
box{18} = [6686 7709; 1413 2122; 235 757]; % 20
box{19} = [547 1554; 2364 3114; 625 1225]; % 38

for i=1:length(box)
    bbox = box{i};
    bboxM4Cropped = round(bsxfun(@plus,bsxfun(@plus,bbox,[1, -2]), ...
             -mag1bbox(:,1))./4);
    bbox = bsxfun(@plus,bsxfun(@plus,4.*bboxM4Cropped,[-1, 2]), ...
             mag1bbox(:,1));
    mask = Util.getBbox(Mask, bboxM4Cropped);
    vessel = Util.getBbox(vessels, bboxM4Cropped);

    BoxList{length(BoxList)+1,1} = {bbox, mask, vessel, root, prefix, length(BoxList)+1};
end
%%
cluster = Cluster.getCluster('-p -400','-tc 30','-l h_vmem=72G','-l s_rt=3:59:00','-l h_rt=4:00:00');
job = Cluster.startJob(@MouseP28.ProcessSingleNucleus, Boxes, 'name','detectNuclei','cluster', cluster) %BoxList
%% Procedure for writing into Knossos Hierarchy, hasn't been done for this dataset
%{
% If two Boxes are overlapping, they have to be processed successively
% by one worker since the parallel writing to webKnossos would generate
% errors.
for i=1:length(BoxList)
    for j=1:length(BoxList)
        matrix(i,j) = NuclearPores.testoverlap(BoxList{i},BoxList{j});
    end
end

vec = 1:length(BoxList);

for i=1:size(matrix,1)
    Corr = NuclearPores.correlationCheck(vec(matrix(i,:)),vec,matrix);
    if i==1
        CorrList{i} = Corr;
    end
    equal = 0;
    for j=1:length(CorrList)
        if isequal(CorrList{j},Corr) == 1
            equal = 1;
        end
    end
    if equal == 0
        CorrList{length(CorrList)+1,1} = NuclearPores.correlationCheck(vec(matrix(i,:)),vec,matrix);
    end
end

for i=1:length(CorrList)
    Boxes{i,1} = {dataset.raw.root,dataset.raw.prefix,{BoxList{CorrList{i}}}};
end

%%
cluster = Cluster.getCluster('-p -400','-tc 50','-l h_vmem=48G','-l s_rt=39:59:00','-l h_rt=40:00:00');
job = Cluster.startJob(@ProcessSingleNucleus, Boxes, 'name','detectNuclei','cluster', cluster)

%%
% (Re)Downsampling KNOSSOS hierachies
dataset.seg.root = strrep(dataset.raw.root,'corrected/color','nuclei/segmentation');
dataset.seg.prefix = strrep(dataset.raw.prefix,'corrected','nuclei');

display('(Re)Downsampling KNOSSOS hierachies for segmentation to include updates');
tic;
thisBBox = [1 1 1; (ceil((mag1bbox_0(:,2)-128)./1024).*1024)']';
createResolutionPyramid(dataset.seg.root, dataset.seg.prefix, thisBBox, strrep(dataset.seg.root, '/1/', ''), true);
toc;
%}
%%
function ProcessSingleNucleus(root, prefix, bboxes)

seg_root = strrep(root,'corrected/color','nuclei/segmentation');
seg_prefix = strrep(prefix,'corrected','nuclei');

gradient_threshold = 0.03;

for i=1:length(bboxes)

    bbox = bboxes{i};
    raw = readKnossosRoi(root, prefix, bbox);
    if min(min(min(raw))) == 0
        boxsize = size(raw)';
        mask = raw==0;
        maskrp = regionprops(permute(~mask, [2 1 3]));
        maskbox = [ceil(maskrp.BoundingBox(1:3))', floor(maskrp.BoundingBox(1:3) + maskrp.BoundingBox(4:6))'];
        raw = Util.getBbox(raw, maskbox);
        bbox(:,1) = bsxfun(@plus, bsxfun(@plus, bbox(:,1), maskbox(:,1)), [-1]);
        bbox(:,2) = bsxfun(@plus, bbox(:,2), -(boxsize - maskbox(:,2)));
    end
    %%
    randoms=randi([110 170],size(raw,1),size(raw,2),size(raw,3));
    darkregion = raw < 60;
    darkregion = imclose(darkregion,strel('disk',4,0));
    for k=1:size(darkregion,3)
        darkregion(:,:,k) = bwareaopen(darkregion(:,:,k), 100);
    end
    darkregion = imdilate(darkregion,strel('disk',10,0));
    
    raw(darkregion) = randoms(darkregion);
    %%
    raw = uint8(smooth3(raw, 'gaussian', 5, 2));
    [Gx, Gy, Gz] = NuclearPores.imgradientxyz(raw);
    gradsum = Gx.^2 + Gy.^2 + Gz.^2;
    clear Gx Gy Gz
    maximum = max(max(max(gradsum)));
    gradsum = gradsum./maximum;
    %%
    nucleus = gradsum > gradient_threshold;
    for k=1:size(nucleus,3)
        nucleus(:,:,k) = bwareaopen(nucleus(:,:,k), 150);
    end
    nucleus = bwareaopen(nucleus, 800);
   
% At this position disturbing compartiments will be erased like
% vessels.
    vessels = readKnossosRoi(strrep(seg_root,'nuclei','vessel'),...
            strrep(seg_prefix,'nuclei','vessel'),bbox, 'uint32');
        
    nucleus = double(nucleus);
    nucleus(logical(vessels))=1;
    %%
    nucleus = imclose(nucleus, NuclearPores.makeSphere(10));
    nucleus = bwareaopen(~nucleus, round(0.01*(sum(sum(sum(nucleus))))));%2e6);
    %%
    darkregion = imdilate(darkregion,NuclearPores.makeSphere(2));
    darkregion = bwareaopen(darkregion,7e5);
    nucleus(darkregion)=1;
    %%
    nucleus = imclose(nucleus,strel('disk',13,0));
    nucleus = padarray(nucleus,[1 1],1);
    nucleus = imfill(nucleus,'holes');
    nucleus = nucleus(2:size(nucleus,1)-1,2:size(nucleus,2)-1,:);

    for k=1:size(nucleus,3)
        nucleus(:,:,k) = ~bwareaopen(~nucleus(:,:,k),3000,4);
    end

    nucleus = bwareaopen(nucleus, round(0.3 * sum(sum(sum(nucleus)))));
    
    if i==1
        nucleus = double(nucleus);
        nucleus(logical(nucleus)) = 2;
        writeKnossosRoi(seg_root, seg_prefix, bbox(:,1)', uint32(nucleus), 'uint32')
    else
        overlap = readKnossosRoi(seg_root, seg_prefix, bbox, 'uint32');
        nucleus = double(nucleus);
        nucleus = nucleus + double(overlap);
        nucleus(logical(nucleus)) = 2;

        writeKnossosRoi(seg_root, seg_prefix, bbox(:,1)', uint32(nucleus), 'uint32')
    end
end
