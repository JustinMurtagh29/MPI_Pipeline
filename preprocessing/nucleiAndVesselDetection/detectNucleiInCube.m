function detectNucleiInCube(datasetMag4, datasetMag4Heur,cubesize,remSmallObjects, offset,newmethod,minmaxArea,visualize,delObjSmallerThan,oldMag4)
if ~exist('newmethod','var') || isempty(newmethod)
    newmethod = true;
end
if ~exist('visualize','var') || isempty(visualize)
    visualize = false;
end
if ~exist('oldMag4','var') || isempty(oldMag4)
    oldMag4 = false;
end
if ~exist('delObjSmallerThan','var') || isempty(delObjSmallerThan)
    delObjSmallerThan = 70000;
end
if ~exist('minmaxArea','var') || isempty(minmaxArea)
    minmaxArea = [3000 45000]; % this is the 2D minimum and maximum area the soma is allowed to have (adjusted for mag4)
end

thisbbox = [offset(:),offset(:)+cubesize-1];
thisbbox = max(1,thisbbox);  % delete negative parts of the bbox
heur = loadRawData(datasetMag4Heur, thisbbox);
raw = loadRawData(datasetMag4, thisbbox);
sizeRaw = size(raw);
nuclei = false(sizeRaw);
if visualize
    figure;%('Position', [3841 1 1920 999]);
    % Create all green RGB image for overlay
end
green = cat(3, zeros([sizeRaw(1) sizeRaw(2)]), ones([sizeRaw(1) sizeRaw(2)]), zeros([sizeRaw(1) sizeRaw(2)]));
for z=1:size(raw,3)
    vessels = heur(:,:,z) == 128;
    % Update for speed at some point
    % Initialization
    thisImage = (raw(:,:,z));
    % Calculate mask (zeros voxel touching border of matrix)
    mask = ~imfill(~imdilate(thisImage == 0, ones([10 10])), 'holes');
    % Detect nuclei
    if newmethod
%         thisImage = imgaussfilt(thisImage, 2);
        thisImage(imclose(thisImage<100,strel('disk',3))) = 0;
        % find regions of similar intensity
        regions = detectMSERFeatures(thisImage, 'ThresholdDelta', 0.5, 'MaxAreaVariation', 0.2, 'RegionAreaRange', minmaxArea);
        % keep only regions which are roughly circular
        div = @(x) x(1)/x(2);
        regions = regions(arrayfun(@(x) div(regions(x).Axes),1:regions.Count) < 2);
        
        if regions.Count > 0
            % keep only regions which have a smooth surrounding boundary (like nuclei)
            % by calculating the convex hull and check how much its area resembles the
            % area of the region
            keepReg = false(regions.Count,1);
            for i=1:regions.Count
                testIm = false(size(thisImage));
                idx = sub2ind(size(thisImage), regions.PixelList{i}(:,2), regions.PixelList{i}(:,1));
                testIm(idx) = 1;
                P = bwboundaries(testIm,'noholes');
                [~,ind] = max(cellfun(@numel,P));
                B = P{ind};
                pa = polyarea(B(:,1),B(:,2));
                [~,ca] = convhull(B(:,1),B(:,2));
                keepReg(i) = ca <= minmaxArea(2) & pa/ca > 0.85 & size(regions.PixelList{i},1)/ca > 0.65;
            end
            regions = regions(keepReg);
            %                 [~,a] = cellfun(@(x) convhull(double(x(:,1)),double(x(:,2))),regions(:).PixelList,'uni',0);
            %                 regions = regions(cellfun(@(x) size(x,1),regions(:).PixelList)./cell2mat(a) > 0.75);
            if regions.Count > 1
                % also check if two regions are for the same nucleus,
                % and if yes, check which of them has a smoother
                % surrounding
                combinations = nchoosek(1:regions.Count,2);
                closeToEachOther = pdist(regions(:).Location) < 25;
                if any(closeToEachOther)
                    regionSets = Graph.findConnectedComponents(combinations(closeToEachOther,:),0,0,regions.Count);
                    
                    [~,a] = cellfun(@(x) convhull(double(x(:,1)),double(x(:,2))),regions(:).PixelList,'uni',0);
                    a = cell2mat(a);
                    [~,ind] = cellfun(@(x) max(a(x)),regionSets);
                    regions = regions(arrayfun(@(x) regionSets{x}(ind(x)),1:numel(regionSets)));
                end
            end
%             figure;
%             imshow(thisImage),hold on
%             plot(regions, 'showPixelList', true, 'showEllipses', true);
        end
        detImage = false(size(thisImage));
        for i=1:regions.Count
            idx = sub2ind(size(detImage), regions.PixelList{i}(:,2), regions.PixelList{i}(:,1));
            detImage(idx) = 1;
        end
        detImage = imopen(detImage,strel('disk',1));
    else
        % Detect apicals for later exclusion
        apicals = medfilt2(thisImage, [5 5]) > 150 & ~vessels;
        apicals = bwareaopen(apicals, 2000);
        apicals = imclose(apicals, strel('disk', 5));
        % detect nuclei with edge detection of nuclei membrane
        detImage = ~edge(thisImage, 'canny');
        detImage = imopen(detImage, strel('disk', 5));
        detImage = detImage & ~apicals;
    end
    detImage = detImage & ~vessels & ~mask;
    detImage = bwareaopen(detImage, 4000, 4);
    detImage = imclose(detImage | mask, strel('disk', 4)) & ~mask;
    detImage = fillHoles(detImage, mask);
    nuclei(:,:,z) = detImage;
    if visualize
        subplot(1,2,1); imshow(raw(:,:,z));
        title(sprintf('Plane %d',z))
        subplot(1,2,2); imshow(raw(:,:,z)); hold on; h = imshow(green); set(h, 'AlphaData', detImage.*0.5);
        title(sprintf('Plane %d',z))
        hold off;
        drawnow;
        pause(0.1)
    end
end

clear raw thisImage apicals
tic;
display('Postprocessing');
% Close to make more constinstent across images
% close with a smaller sphere than opening
if oldMag4
    % in mag 4 this creates a sphere of ~250 nm radius
    [x,y,z] = meshgrid(-4:4,-4:4,-2:2);
    se1 = (x/4).^2 + (y/4).^2 + (z/2).^2 <= 1;
    % in mag 4 this creates a sphere of ~450 nm radius
    [x,y,z] = meshgrid(-10:10,-10:10,-4:4);
    se2 = (x/10).^2 + (y/10).^2 + (z/4).^2 <= 1;
else
    % in mag 4 this creates a sphere of ~250 nm radius
    [x,y,z] = meshgrid(-4:4,-4:4,-3:3);
    se1 = (x/4).^2 + (y/4).^2 + (z/2).^2 <= 1;
    % in mag 4 this creates a sphere of ~450 nm radius
    [x,y,z] = meshgrid(-10:10,-10:10,-7:7);
    se2 = (x/10).^2 + (y/10).^2 + (z/7).^2 <= 1;
end
nuclei = imclose(nuclei, se1); % merge soma parts where a few z planes are missing or some hole occurred
nuclei = imopen(nuclei, se2); % divide merged somata (e.g. when label leaked into cytosol and then across soma membranes)
if remSmallObjects
    % Remove very small objects (smaller than 100*70*10 voxel
    nuclei = bwareaopen(nuclei, delObjSmallerThan);
end
if 0
    figure;
    for z=1:size(raw,3)
        imshow(raw(:,:,z)); hold on; h = imshow(green); set(h, 'AlphaData', nuclei(:,:,z).*0.5);
        title(sprintf('Plane %d',z))
        hold off;
        drawnow;
        pause(0.3)
    end
end
    
heur(nuclei) = 255;

saveRawData(datasetMag4Heur, max(1,offset),heur);  % save and again delete negative parts

toc;
