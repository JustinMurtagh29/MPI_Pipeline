addpath(genpath('/gaba/u/alik/code/pipeline'));

root = '/mnt/disk-two/V1/bvn/4/';
prefix = 'V1_D_JK_18052017_mag4';
myDataset=KnossosDataset(root)
bbox=myDataset.getBbox();
raw = readKnossosRoi(root, prefix, bbox);
save('raw','raw','-v7.3');
output = '/mnt/disk-two/V1/bvn/';
%% Detect vessels with visualization flag set to false
vessels = AKdetectVessels(raw, true,output); % 7000s on 07x2 (mostly postprocessing)
save(fullfile(output,'vesselsFinal'),'vessels','-v7.3');
makeSegMovie(vessels, raw, fullfile(output,'BVmovieFinal.avi'), 1);
%% Give each detected vessel CC its own unique ID
vesselL = labelmatrix(bwconncomp(vessels));
rp = regionprops(vesselL);
for i =[2,4,6,7]
vessels(vesselL == i) = 0;
end
%% Visualize result with raw data for weird looking isosurfaces
% [9 68 77 96 111 117]
for i=7
    for z=ceil(rp(i).BoundingBox(3)):ceil((rp(i).BoundingBox(3)+rp(i).BoundingBox(6)))
        x = floor(rp(i).BoundingBox(2)):ceil(rp(i).BoundingBox(2)+rp(i).BoundingBox(5));
        y = floor(rp(i).BoundingBox(1)):ceil(rp(i).BoundingBox(1)+rp(i).BoundingBox(4));
        imshow(raw(x,y,z)); hold on; 
        temp = label2rgb(vessels(x,y,z), [0 1 0; 1 0 0; 0 0 1], [1 1 1]);
        himage = imshow(temp);
        set(himage, 'AlphaData', 0.3);
        title(num2str(i));
        pause(1);
        hold off;
    end
end
%% Visualize all CC individually to judge/search errors
a = vesselL(1:6:end,1:6:end,1:end);
for i=7
    close all; figure;
    isosurface(a == i, .5);
    drawnow;
    
    camorbit(0,180); camlight;
   
    camorbit(90,0); camlight;

end


vessels(both) = true;
nuclei(both) = false;
%% give each CC a unique label
nucleiL = labelmatrix(bwconncomp(nucleiPost));
rp = regionprops(nucleiL);
%% nuclei raw/nuclei flag visualization
%[62,72,122,128,141,160]
for i=1:size(rp,1)
   
    for z=ceil(rp(i).BoundingBox(3)):2:floor((rp(i).BoundingBox(3)+rp(i).BoundingBox(6)))
        x = ceil(rp(i).BoundingBox(2)):floor(rp(i).BoundingBox(2)+rp(i).BoundingBox(5));
        y = ceil(rp(i).BoundingBox(1)):floor(rp(i).BoundingBox(1)+rp(i).BoundingBox(4));
        imshow(raw(x,y,z)); hold on; 
        temp = label2rgb(nucleiPost(x,y,z), [0 1 0; 1 0 0; 0 0 1], [1 1 1]);
        himage = imshow(temp);
        set(himage, 'AlphaData', 0.3);
        title(num2str(i));
        pause(.1);
        hold off;
    end
end

%% deleting wrongconnected components
for i=[36 45 123 135 142]
nucleiPost(nucleiL == i) = 0;
end
    %% 62 some edges added because of the yz detection deleted here
i = 81;
x = floor(rp(i).BoundingBox(2)):floor(rp(i).BoundingBox(2)+rp(i).BoundingBox(5));
y = floor(rp(i).BoundingBox(1)):floor(rp(i).BoundingBox(1)+rp(i).BoundingBox(4));
z = floor(rp(i).BoundingBox(3)):floor(rp(i).BoundingBox(3)+rp(i).BoundingBox(6));
isosurface(nucleiPost(x,y,z), .5);

roi = nucleiPost(x,y,z);
roi(:,:,83:86)=0;
roi=bwareaopen(roi,30000);
figure;isosurface(roi,.5)
% Fix in original data
nucleiPost(x,y,z) = roi;
%save(fullfile('nuclei80fixed'),'nucleiPost','-v7.3');
%% 
i = 81;
x = floor(rp(i).BoundingBox(2)):floor(rp(i).BoundingBox(2)+rp(i).BoundingBox(5));
y = floor(rp(i).BoundingBox(1)):floor(rp(i).BoundingBox(1)+rp(i).BoundingBox(4));
z = floor(rp(i).BoundingBox(3)):floor(rp(i).BoundingBox(3)+rp(i).BoundingBox(6));
figure();isosurface(nucleiPost(x,y,z), .5);
daspect([1 1 1])


% Test which size of structuring element we need to make it decompose into
% two parts
roi = nucleiPost(x,y,z);
distTraf = -bwdist(~roi);
minima = imextendedmin(distTraf, 1);
distTraf = imimposemin(distTraf, ~roi | minima);
fixed = watershed(distTraf);

% Remove background object
fixed(fixed == 1) = 0;
% Seperate objects
roiNew = imclose(fixed == 4 , ones(3,3,3)) | imclose(fixed == 3 |fixed == 2 , ones(3,3,3));

for i=1:size(fixed,3)
    imagesc(fixed(:,:,i));
    pause(.5);
end

for i=1:size(roiNew,3)
    imshow(roiNew(:,:,i));
    pause(.5);
end

% Fix in original data
%nucleiPost(x,y,z) = roiNew;
%% Save
makeSegMovie(nuclei, raw, fullfile(output,'nucleiMag4.avi'), 1);
save(fullfile(output,'nucleiMag4Finished'),'nucleiPost','-v7.3');

%% Joint isosurface

figure;
temp = nucleiPost(1:6:end,1:6:end,:);
n = patch(isosurface(temp, .5)); hold on;
temp = vessels(1:6:end,1:6:end,:);
v = patch(isosurface(temp, .5));

n.FaceColor = 'red';
v.FaceColor = 'blue';
n.EdgeColor = 'none';
v.EdgeColor = 'none';
daspect([28 28 11.24*6]);

view(3); axis tight
camlight 
lighting gouraud
saveas(gcf, '/home/mberning/Desktop/both.png');

%%
n = regionprops(nucleiPost);
v = regionprops(vessels);

figure;
hist([n(:).Area].*prod([28 28 11.24].*4)/1000^3, 50);
title('Nucleus size [microns^3]');

%%
rawNuclei = raw(nuclei);
regionprops(rawNuclei, 'all');

%% Give each nucleus its own unique ID
nucleiL = labelmatrix(bwconncomp(nucleiPost));
rp = regionprops(nucleiL);

%% Visualize all CC individually to judge/search errors
a = nucleiL(1:6:end,1:6:end,1:end);
for i=1:max(a(:))
    close all; figure;
    isosurface(a == i, .5);
    drawnow;
    saveas(gcf,['/home/mberning/Desktop/nuclei/' num2str(i, '%.3i') '_1.png']);
    camorbit(0,180); camlight;
    saveas(gcf,['/home/mberning/Desktop/nuclei/' num2str(i, '%.3i') '_2.png']);
    camorbit(90,0); camlight;
    saveas(gcf,['/home/mberning/Desktop/nuclei/' num2str(i, '%.3i') '_3.png']);
end



%% SAME POSTPROCESSING STARTS HERE, remove errors in detection

%% 9 is cell with two nuclei (glia) and one merged endothelial nucleus
i = 127;
x = floor(rp(i).BoundingBox(2)):ceil(rp(i).BoundingBox(2)+rp(i).BoundingBox(5));
y = floor(rp(i).BoundingBox(1)):ceil(rp(i).BoundingBox(1)+rp(i).BoundingBox(4));
z = floor(rp(i).BoundingBox(3)):ceil(rp(i).BoundingBox(3)+rp(i).BoundingBox(6));
isosurface(nucleiPost(x,y,z), .5);

% Test which size of structuring element we need to make it decompose into
% two parts
roi = nucleiPost(x,y,z);
distTraf = -bwdist(~roi);
minima = imextendedmin(distTraf, 1);
distTraf = imimposemin(distTraf, ~roi | minima);
fixed = watershed(distTraf);

for i=1:3:size(fixed,2)
    imagesc(squeeze(fixed(:,i,:)));
    pause(.5);
end

% Remove background object
fixed(fixed == 1) = 0;
% Seperate objects
%roiNew = imclose(fixed == 2 | fixed == 3,ones(3,3,3)) | imclose(fixed == 4 | fixed == 5, ones(3,3,3));
roiNew = imclose(fixed == 2 ,ones(3,3,3)) | imclose(fixed == 3|fixed == 4| fixed == 5, ones(3,3,3));
% Fix in original data
nucleiPost(x,y,z) = roiNew;

for i=1:3:size(roiNew,2)
    imshow(squeeze(roiNew(:,i,:)));
    pause(.5);
end


%% 68 apical -> delete
nuclei(nucleiL == 68) = 0;

%% 111 merged blood vessel part
i = 111;
x = floor(rp(i).BoundingBox(2)):ceil(rp(i).BoundingBox(2)+rp(i).BoundingBox(5));
y = floor(rp(i).BoundingBox(1)):ceil(rp(i).BoundingBox(1)+rp(i).BoundingBox(4));
z = floor(rp(i).BoundingBox(3)):ceil(rp(i).BoundingBox(3)+rp(i).BoundingBox(6));

roi = nuclei(x,y,z);
distTraf = -bwdist(~roi);
minima = imextendedmin(distTraf, .7);
distTraf = imimposemin(distTraf, ~roi | minima);
fixed = watershed(distTraf);

for i=1:size(fixed,3)
    imagesc(fixed(:,:,i));
    title(num2str(i));
    pause(1);
end

fixed(fixed == 1) = 0;
fixed(fixed == 2) = 0;
roiNew = imclose(fixed == 3, ones(3,3,3));
nuclei(x,y,z) = roiNew;

%% 117 endothelial nucleus, but some apical attachted
i = 117;
x = floor(rp(i).BoundingBox(2)):ceil(rp(i).BoundingBox(2)+rp(i).BoundingBox(5));
y = floor(rp(i).BoundingBox(1)):ceil(rp(i).BoundingBox(1)+rp(i).BoundingBox(4));
z = floor(rp(i).BoundingBox(3)):ceil(rp(i).BoundingBox(3)+rp(i).BoundingBox(6));

roi = nuclei(x,y,z);
distTraf = -bwdist(~roi);
minima = imextendedmin(distTraf, 1);
distTraf = imimposemin(distTraf, ~roi | minima);
fixed = watershed(distTraf);

for i=1:size(fixed,3)
    imagesc(fixed(:,:,i));
    pause(.5);
end

fixed(fixed == 1) = 0;
fixed(fixed == 3) = 0;
roiNew = imclose(fixed == 2, ones(3,3,3));
nuclei(x,y,z) = roiNew;
%%
nucleiL = labelmatrix(bwconncomp(nucleiPost));
rp = regionprops(nucleiL);

save('nucleiL','nucleiL','rp','-v7.3')
%%
writeKnossosRoi('/run/media/mberning/localStorage/07x2_mag4_vessels/', 'vessels', [1 1 1], vessels);
writeKnossosRoi('/run/media/mberning/localStorage/07x2_mag4_nuclei/', 'nuclei', [1 1 1], nucleiL);
