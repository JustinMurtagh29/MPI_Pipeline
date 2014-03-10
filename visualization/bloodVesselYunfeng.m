path = 'C:\Users\mberning\Desktop\yunfeng\';
files = dir([path '*.jpg']);
% Read images
for i=1:length(files)
    img{i} = imread([path files(i).name]);  
end
% Choose ROI
img{1} = img{1}(51:end-50,91:end-10,:);
img{2} = img{2}(11:end-90,51:end-50,:);
img{3} = img{3}(61:end-40,110:end-41,:);
img{4} = img{4}(52:end-51,76:end-75,:);
img{5} = img{5}(60:end-159,12:end-112,:);
img{6} = img{6}(12:end-91,36:end-115,:);
img{7} = img{7}(17:end-116,13:end-113,:);
img{8} = img{8}(68:end-107,11:end-90,:);
% Convert to grayscale and equalize all histograms
for i=1:length(files)
    img{i} = rgb2gray(img{i});
    img{i} = histeq(img{i}, 256);
end
% Different hemisphere, mirror
for i=5:8
    img{i} = img{i}(1:end,end:-1:1);
end

%% Look at pictures
close all;
for i=1:length(files)
    figure; imshow(img{i});
end

%%
barrelPos = [216 151; 204 161; 279 213; 183 165; 264 162; 215 142; 192 165; 282 169];
barrelPos(:,[1 2]) = barrelPos(:,[2 1]);
meanBarrel = round(mean(barrelPos,1));
center = [1000 1000];
sumImage = zeros(2000,2000);
alignedImage = zeros(2000,2000,8, 'uint8');
left = 0;
bottom = 0;
right = 1000;
top = 1000;
for i=1:length(img)
    x = bsxfun(@plus, center-floor(size(img{i},1)/2):center+floor(size(img{i},1)/2)-1, meanBarrel(1)-barrelPos(i,1));
    y = bsxfun(@plus, center-floor(size(img{i},2)/2):center+floor(size(img{i},2)/2)-1, meanBarrel(2)-barrelPos(i,2));
    sumImage(x,y) = sumImage(x,y) +  single(img{i});
    left = max(left,min(x));
    bottom = max(bottom,min(y));
    right = min(right,max(x));
    top = min(top,max(y));
    alignedImage(x,y,i) = img{i};
end
implay(alignedImage(left:right,bottom:top,:));
% imagesc(sumImage(left:right,bottom:top));
% colormap('gray');