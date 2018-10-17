function nucleiPost = AKdetectNuclei( raw, vessels, visualize )
% Pass raw data array (tested with 07x2 mag4) and detect somata in dataset,
% pass visualize flag in order to tune parameter on new dataset, goal
% should be to detect Somata with close to 0 FP rate, some FN are no problem, due to
% closing step in postprocessing, could be changed by adding opening

% By default: visualization off
if nargin < 3
    visualize = false;
end

sizeRaw = size(raw);
nuclei = false(size(raw));

if visualize
    figure('Position', [3841 1 1920 999]);
    % Create all green RGB image for overlay
    green = cat(3, zeros([sizeRaw(1) sizeRaw(2)]), ones([sizeRaw(1) sizeRaw(2)]), zeros([sizeRaw(1) sizeRaw(2)]));
end

tic;
for z=1:sizeRaw(3)
    % Update for speed at some point
    display(['Processing slice: ' num2str(z)]);
    % Initialization
    thisImage = single(raw(:,:,z));
    theseVessels = vessels(:,:,z);
    % Calculate mask (zeros voxel touching border of matrix)
    mask = ~imfill(~imdilate(thisImage == 0, ones([10 10])), 'holes');
    % Detect apicals for later exclusion
    apicals = medfilt2(thisImage, [5 5]) > 150 & ~theseVessels;
    apicals = bwareaopen(apicals, 2000);
    apicals = imclose(apicals, strel('disk', 5));
    % Detect nuclei
    thisImage = ~edge(thisImage, 'canny');
    thisImage = imopen(thisImage, strel('disk', 5));
    theseNuclei = thisImage & ~vessels(:,:,z) & ~mask & ~apicals;
    theseNuclei = bwareaopen(theseNuclei, 2000, 4);
    theseNuclei = imclose(theseNuclei | mask, strel('disk', 3)) & ~mask;
    theseNuclei = fillNucleiHoles(theseNuclei, mask);
    % Save for output
    nuclei(:,:,z) = theseNuclei;
    if visualize 
        subplot(1,2,1); imshow(raw(:,:,z));
        hold off;
        subplot(1,2,2); imshow(raw(:,:,z)); hold on; h = imshow(green); set(h, 'AlphaData', theseNuclei.*0.2);
        hold off;
        drawnow;
    end
end
toc;

tic;
display('Postprocessing');
clear raw vessels
% Close to make more constinstent across images
[x,y,z] = meshgrid(-8:8,-8:8,-3:3);
se = (x/8).^2 + (y/8).^2 + (z/3).^2 <= 1;
nucleiPost = imclose(nuclei, se);
nucleiPost = imopen(nucleiPost, se);
% Remove very small objects (smaller than 100*30*10 voxel
nucleiPost = bwareaopen(nucleiPost, 30000);
toc;

end

function nuclei = fillNucleiHoles(nuclei, mask)
    % Dataset specific changes needed here, this function takes mask and
    % drill holes into 'outer hull' to be able to use imfill, change drill
    % location according to dataset :)
    vesselMasked = or(nuclei,mask);
    vesselMasked(1:400, 800) = 0;
    nuclei = imfill(vesselMasked, 'holes') & ~mask;
end