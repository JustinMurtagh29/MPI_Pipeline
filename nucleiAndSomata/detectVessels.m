function vesselsPost = detectVessels( raw, visualize )
% Pass raw data array (tested with 07x2 mag4) and detect vessel in dataset,
% pass visualize flag in order to tune parameter on new dataset, goal
% should be to detect vessels with 0 FP rate, some FN are no problem, due to
% closing step in postprocessing

% By default: visualization off
if nargin == 1
    visualize = false;
end

sizeRaw = size(raw);
vessels = false(size(raw));

if visualize
    figure('Position', [3841 1 1920 999]);
    % Create all green RGB image for overlay
    green = cat(3, zeros([sizeRaw(1) sizeRaw(2)]), ones([sizeRaw(1) sizeRaw(2)]), zeros([sizeRaw(1) sizeRaw(2)]));
end

tic;
for z=1:sizeRaw(3)
    display(['Processing slice: ' num2str(z)]);
    thisImage = raw(:,:,z);
    theseVessels = false(size(thisImage));
    % Calculate mask (zeros voxel touching border of matrix)
    mask = ~imfill(~imdilate(thisImage == 0, ones([10 10])), 'holes');
    % Median Filter on image & keep only high pixel values
    thisImage = medfilt2(thisImage, [3 3]);
    thisImage(~mask & thisImage < 150) = 150;
    % These are only detection parameter that may need tuning on new dataset or
    % different magnification level (in addition to parameter in line above)
    regions = detectMSERFeatures(thisImage, 'ThresholdDelta', 2, 'MaxAreaVariation', 0.5, 'RegionAreaRange', [1000 100000]);
    % Set regions to true that were detected by MSER
    for i=1:regions.Count
        idx = sub2ind(size(thisImage), regions.PixelList{i}(:,2), regions.PixelList{i}(:,1));
        theseVessels(idx) = 1;
    end
    % Fill holes in detection
	theseVessels = fillVesselHoles(theseVessels, mask);
    % Save for output
    vessels(:,:,z) = theseVessels;
    if visualize
        subplot(1,2,1); imshow(raw(:,:,z));
        subplot(1,2,2); imshow(raw(:,:,z)); hold on; h = imshow(green); set(h, 'AlphaData', theseVessels.*0.5);
        hold off;
        drawnow;
    end
end
toc;

tic;
display('Postprocessing');
% Smooth and connect surface (by closing with spherical (in nm) strucuturing element)
[x,y,z] = meshgrid(-12:12,-12:12,-5:5);
se = (x/12).^2 + (y/12).^2 + (z/5).^2 < 1;
vesselsPost = imclose(vessels, se);
% Remove very small objects (smaller than 100*10*10 voxel
vesselsPost = bwareaopen(vesselsPost, 10000);
toc;

end

function vessel = fillVesselHoles(vessel, mask)
    % Dataset specific changes needed here, this function takes mask and
    % drill holes into 'outer hull' to be able to use imfill, change drill
    % location according to dataset :)
    vesselMasked = or(vessel,mask);
    vesselMasked(1:400, 200) = 0;
    vessel = imfill(vesselMasked, 'holes') & ~mask;
end