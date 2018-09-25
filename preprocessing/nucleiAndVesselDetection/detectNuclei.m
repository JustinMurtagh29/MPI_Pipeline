function nuclei = detectNuclei( raw, vessels, visualize,remSmallObjects )
% Pass raw data array (tested with 07x2 mag4) and detect somata in dataset,
% pass visualize flag in order to tune parameter on new dataset, goal
% should be to detect Somata with close to 0 FP rate, some FN are no problem, due to
% closing step in postprocessing, could be changed by adding opening

% By default: visualization off
if nargin < 3 || isempty(visualize)
    visualize = false;
end
if nargin < 4 || isempty(remSmallObjects)
    remSmallObjects = true;
end

sizeRaw = size(raw);
nuclei = false(size(raw));

if visualize
    figure('Position', [3841 1 1920 999]);
end
% Create all green RGB image for overlay
green = cat(3, zeros([sizeRaw(1) sizeRaw(2)]), ones([sizeRaw(1) sizeRaw(2)]), zeros([sizeRaw(1) sizeRaw(2)]));
tic;
parfor z=1:sizeRaw(3)  % parfor takes too much memory..
    % Update for speed at some point
    display(['Processing slice: ' num2str(z)]);
    % Initialization
    thisImage = single(raw(:,:,z));
    % Calculate mask (zeros voxel touching border of matrix)
    mask = ~imfill(~imdilate(thisImage == 0, ones([10 10])), 'holes');
    % Detect apicals for later exclusion
    apicals = medfilt2(thisImage, [5 5]) > 150 & ~vessels(:,:,z);
    apicals = bwareaopen(apicals, 2000);
    apicals = imclose(apicals, strel('disk', 5));
    % Detect nuclei
    thisImage = ~edge(thisImage, 'canny');
    thisImage = imopen(thisImage, strel('disk', 5));
    thisImage = thisImage & ~vessels(:,:,z) & ~mask & ~apicals;
    thisImage = bwareaopen(thisImage, 2000, 4);
    thisImage = imclose(thisImage | mask, strel('disk', 3)) & ~mask;
    thisImage = fillHoles(thisImage, mask);
    % Save for output
    nuclei(:,:,z) = thisImage;
    if visualize 
        subplot(1,2,1); imshow(raw(:,:,z));
        hold off;
        subplot(1,2,2); imshow(raw(:,:,z)); hold on; h = imshow(green); set(h, 'AlphaData', thisImage.*0.2);
        hold off;
        drawnow;
    end
end
toc;
clear raw vessels thisImage apicals
tic;
display('Postprocessing');
% Close to make more constinstent across images
[x,y,z] = meshgrid(-8:8,-8:8,-3:3);
se = (x/8).^2 + (y/8).^2 + (z/3).^2 <= 1;
nuclei = imclose(nuclei, se);
nuclei = imopen(nuclei, se);
if remSmallObjects
    % Remove very small objects (smaller than 100*30*10 voxel
    nuclei = bwareaopen(nuclei, 30000);
end
toc;

end