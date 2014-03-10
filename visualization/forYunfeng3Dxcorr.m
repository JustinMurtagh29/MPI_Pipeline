% For Yunfeng. 3D cross correlation
load('I:\CortexConnectomics\Yunfeng\bv_2p.mat');
load('I:\CortexConnectomics\Yunfeng\bv_em.mat');
% Create intiutive names
em = bw;
lm = data;
% Voxel Sizes estimated
voxelSizeLm = [1.3 1.3];
voxelSizeEm = [1/4.9 1/4.9];
% Downsample em data
downscale = ceil(voxelSizeLm./voxelSizeEm);
emDS = nlfilter3(single(em(1:end-mod(size(em,1),downscale(1)),1:end-mod(size(em,1),downscale(2)),:)), @max, [downscale 1]);
emDS = padarray(emDS, [100 100 0], 0);
emDS(emDS == 0) = -1;
display('Data loaded and preprocessed');

%% Rotate data using NN interpolation
rotations = 0:45:315;
dataSize = size(emDS);
border = ([(dataSize(1:2) - [1 1])/2 0]);
for i=1:length(rotations)
    orth1New = rotate([1 0 0], [0 0 1], rotations(i));
    orth2New = rotate([0 1 0], [0 0 1], rotations(i));
    % Define grid for original values (emDS)
    [X,Y,Z] = meshgrid(1-border(1):dataSize(1)-border(1),1-border(2):dataSize(2)-border(2),1:dataSize(3));
    % Construct rotation matrix
    A = [orth1New; orth2New; [0 0 1]];
    % Define grid for interpolating values (involves rotation)
    [XI, YI, ZI] = meshgrid(1-border(1):dataSize(1)-border(1),1-border(2):dataSize(2)-border(2),1:dataSize(3)); %new grid in new coordinate system
    coords = [XI(:) YI(:) ZI(:)];
    coords = (inv(A)*coords')'; % new grid in old coordinate system
    XI = reshape(coords(:,1),size(XI));
    YI = reshape(coords(:,2),size(YI));
    ZI = reshape(coords(:,3),size(ZI));
    % Do the interpolation
    emRotated{i} = interp3(X, Y, Z, emDS, XI, YI, ZI, 'nearest');
    emRotated{i}(isnan(emRotated{i})) = -1;
end

% %% Look at it
% for i=1:length(rotations)
%     isosurface(emRotated{i}, 0);
%     drawnow;
%     pause(.5);
% end
% close all;

%% Look at it
figure('Position', [1 1 1600 1124]);
rot = 1;
caxis([0 max(emRotated{rot}(:))]);
for i=1:size(emRotated{rot},3)
    imagesc(emRotated{rot}(:,:,i));
    colorbar
    pause(.1);
end
close all;

%% Cross-Correlate
matlabpool 4;
lm = single(lm);
lm(lm == 0) = -1;
parfor rot=1:length(rotations)
    tic
    crossCorrelation{rot} = convn(lm, emRotated{rot}(end:-1:1,end:-1:1,end:-1:1), 'full');
    t = toc;
    display([num2str(rot) '. cross correlation done: ' num2str(t, '%6.0f') ' seconds.']);
end
matlabpool close;
save('D:\Yunfeng\firstTryCrossCorrelationForYH.mat', '-v7.3');

%% Find maximum
unvalidBorder = min(size(lm),size(emDS)) + [102 102 0]; 
for i=1:length(rotations)
    crossTemp = crossCorrelation{i}(1+unvalidBorder(1):end-unvalidBorder(1),1+unvalidBorder(2):end-unvalidBorder(2),1+unvalidBorder(3):end-unvalidBorder(3));
    [v(i), idx(i)] = max(crossTemp(:));
    [x(i) y(i) z(i)] = ind2sub(size(crossTemp), idx(i)); % Position of maximum in each image
end
[vR, idxR] = max(v);
xTest = x(idxR) + unvalidBorder(1);
yTest = y(idxR) + unvalidBorder(2);
zTest = z(idxR) + unvalidBorder(3);

%% Look at it
figure('Position', [1 1 1600 1124]);
rot = idxR;
caxis([min(min(min(crossCorrelation{rot}(1+unvalidBorder(1):end-unvalidBorder(1),1+unvalidBorder(2):end-unvalidBorder(2),1+unvalidBorder(3):end-unvalidBorder(3))))) ...
    max(max(max(crossCorrelation{rot}())))]);
for i=1+unvalidBorder(3):size(crossCorrelation{rot},3)-unvalidBorder(3)
    imagesc(crossCorrelation{rot}(1+unvalidBorder(1):end-unvalidBorder(1),1+unvalidBorder(2):end-unvalidBorder(2),i));
    axis equal; axis off;
    pause(.1);
end
close all;

%% Calculate isosurfaces
issf = isosurface(smooth3(single(lm), 'gaussian', [5 5 5], 1), 0);
issfEM = isosurface(smooth3(emRotated{idxR}, 'gaussian', [5 5 5], 1), 0);
xS = yTest - size(emDS,2);
yS = xTest - size(emDS,1);
zS = zTest - size(emDS,3);
issfEM.vertices = bsxfun(@plus, issfEM.vertices, [yS xS zS]);

%% Isosurfaces of both (at position of maximum cross correlation
figure('Position', [1 1 1600 1124], 'Renderer', 'OpenGL');
hpatch = patch(issf);
set(hpatch,'FaceColor', 'g', 'EdgeColor', 'none');
set(hpatch,'AmbientStrength', .2, 'SpecularStrength',.6, 'DiffuseStrength',.6);
hpatch = patch(issfEM);

set(hpatch,'FaceColor', 'r', 'EdgeColor', 'none');
set(hpatch,'AmbientStrength', .2, 'SpecularStrength',.6, 'DiffuseStrength',.6);
axis equal; axis off;
lighting gouraud;
hlight = camlight('headlight');
daspect([1.3 1.3 2]);
view(3);
