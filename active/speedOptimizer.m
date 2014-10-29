profile on;
for i=1000:1001
    tic;
    object = seg == edges2(i,1);
    otherObject = seg == edges2(i,2);
    border = bwdist(object, 'euclidean') <= 1 & bwdist(otherObject, 'euclidean') <= 1;
    toc;
end
profile viewer
profile off

%%
tic;
temp2 = aff(border);
toc;

affFeatures = [aff calcFeatures(aff)]

%%
tic;
properties = {'BoundingBox' 'Centroid' 'PixelList', 'Area'};
propObject = regionprops(object, properties);
pc = princomp(propObject.PixelList, 'econ');
directionObject = pc(:,1);
toc;