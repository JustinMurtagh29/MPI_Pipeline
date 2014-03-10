objProps = regionprops(seg, 'Centroid');

%%
figure;
hold on;
for i=1:length(edges)
    X = [objProps(edges(i,1)).Centroid(1) objProps(edges(i,2)).Centroid(1)];
    Y = [objProps(edges(i,1)).Centroid(2) objProps(edges(i,2)).Centroid(2)];
    Z = [objProps(edges(i,1)).Centroid(3) objProps(edges(i,2)).Centroid(3)];
    plot3(X,Y,Z, 'LineWidth', 1, 'Color', [1-p(i) p(i) 0]); 
end
%% to Amira
pos = cat(1,objProps(:).Centroid);
pos = pos .* repmat([11.28 11.28 28], size(pos,1), 1);
convertKnossosNmlToHoc3(edges, pos, p, 0, 1, 'C:\Users\mberning\Desktop\allEdges');
convertKnossosNmlToHoc3(edges, pos, p, 0, 0.15, 'C:\Users\mberning\Desktop\rejected');
convertKnossosNmlToHoc3(edges, pos, p, 0.15, 0.95, 'C:\Users\mberning\Desktop\querry');
convertKnossosNmlToHoc3(edges, pos, p, 0.95, 1, 'C:\Users\mberning\Desktop\accepted');

%% Generate isosurfaces

for i=1:max(seg(:))
    temp = seg == i;
    xIdx = any(any(temp,2),3);
    x = [find(xIdx, 1, 'first') find(xIdx, 1, 'last')];
    yIdx = any(any(temp,1),3);
    y = [find(yIdx, 1, 'first') find(yIdx, 1, 'last')];
    zIdx = any(any(temp,1),2);
    z = [find(zIdx, 1, 'first') find(zIdx, 1, 'last')];
    temp = temp(x(1):x(2),y(1):y(2),z(1):z(2));
    if length(size(temp)) == 3
        issfs{i} = isosurface(temp, .5);
        issfs{i}.vertices = bsxfun(@plus, issfs{i}.vertices, [y(1)-1 x(1)-1 z(1)-1]);
    end
    display(num2str(i,'%.2i'));
end

%% bulls***
raw = readKnossosRoi('Z:\CortexConnectomics\shared\cortex\2012-09-28_ex145_07x2\8\', '2012-09-28_ex145_07x2_mag8', [0 1250; 0 1250; 0 500]);
temp = raw ~= 0;

xIdx = sum(sum(~temp,2),3);
yIdx = sum(sum(~temp,1),3);
zIdx = sum(sum(~temp,1),2);

%% stacks aus levelcreator zu video
img = imread('C:\Users\mberning\Downloads\51b32345e4b0be083364a97a\stack.png');

%%
load('C:\Users\mberning\Desktop\issfsForFS.mat');
load('C:\sync\forFS\labelcolorsToMatlab.mat');

idx = zeros(length(issfs),1);
for i=1:length(issfs)
        if isempty(issfs{i}) || isempty(issfs{i}.vertices)
                idx(i) = 1;
        end
end

cm(:,4) = [];
cm(length(issfs)+1:end,:) = [];
issfs(find(idx)) = [];
cm(find(idx),:) = [];

for i=1:length(issfs)
    nrVert(i) = length(issfs{i}.vertices);
end

issfs(nrVert < 20000 | nrVert > 50000) = [];
cm(nrVert < 20000 | nrVert > 50000,:) = [];

for i=1:length(issfs)
    issfs{i}.vertices = issfs{i}.vertices .* repmat([11.28 11.28 28], size(issfs{i}.vertices,1),1);
end

for i=1:length(issfs)
    issfs{i}.vertices(:,[1 2]) = issfs{i}.vertices(:,[2 1]);
end

exportSurfaceToAmira(issfs, 'C:/sync/forFS/smallObjects', cm);

%%
temp = zeros(size(seg), 'uint8');
temp(segNew == 18797) = 1;
temp(segNew == 20416) = 2;
xIdx = any(any(temp,2),3);
x = [find(xIdx, 1, 'first') find(xIdx, 1, 'last')];
yIdx = any(any(temp,1),3);
y = [find(yIdx, 1, 'first') find(yIdx, 1, 'last')];
zIdx = any(any(temp,1),2);
z = [find(zIdx, 1, 'first') find(zIdx, 1, 'last')];

temp = temp(x(1):x(2),y(1):y(2),z(1):z(2));
figure;
isosurface(smooth3(temp == 1, 'gaussian',5,3));
hold on;
isosurface(smooth3(temp == 2, 'gaussian',5,3));
border = imdilate(temp == 1, ones(3,3,3)) & imdilate(temp == 2, ones(3,3,3));
isosurface(border);

%%
a = patch([31280 31280 35000 35000],[0 1 1 0],'g', 'facecolor',[0 1 0], 'edgecolor',[0 1 0], 'facealpha' ,0.5);
b = patch([24080 24080 31280 31280],[0 1 1 0],'b', 'facecolor',[0 0 1], 'edgecolor',[1 1 1], 'facealpha' ,0.5);
c = patch([0 0 24080 24080],[0 1 1 0],'r', 'facecolor',[1 0 0], 'edgecolor',[1 0 0], 'facealpha' ,0.5);

