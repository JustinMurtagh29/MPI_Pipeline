function vessels = detectVesselsApicalsMask(raw, visualize)%, tunnelCoord, threshold)
sizeRaw = size(raw);
vessels = false(size(raw));
threshold = 200;
if visualize
    figure%('Position', [3841 1 1920 999]);
    % Create all green RGB image for overlay
    green = cat(3, zeros([sizeRaw(1) sizeRaw(2)]), ones([sizeRaw(1) sizeRaw(2)]), zeros([sizeRaw(1) sizeRaw(2)]));
end

tic;
for z=1:sizeRaw(3)
    display(['Processing slice: ' num2str(z)]);
    thisImage = raw(:,:,z);                 
%     theseVessels = false(size(thisImage));
    % Calculate mask (zeros voxel touching border of matrix)
    mask = ~imfill(~imdilate(thisImage == 0, ones([10 10])), 'holes');
    % Median Filter on image & keep only high pixel values
    thisImage = medfilt2(thisImage, [3 3]);
    thisImage = ~mask & thisImage > threshold;%(~mask & thisImage < threshold) = threshold;
    % These are only detection parameter that may need tuning on new dataset or
    % different magnification level (in addition to parameter in line above)
%     regions = detectMSERFeatures(thisImage, 'ThresholdDelta', 2, 'MaxAreaVariation', 0.5, 'RegionAreaRange', [1000 100000]);
    % Set regions to true that were detected by MSER
%     for i=1:regions.Count
%         idx = sub2ind(size(thisImage), regions.PixelList{i}(:,2), regions.PixelList{i}(:,1));
%         theseVessels(idx) = 1;
%     end
    theseVessels = bwareaopen(thisImage, 500,8);
    theseVessels = imclose(theseVessels, strel('disk', 20));
    % Fill holes in detection
    tunnelCoord = 1350;
    theseVessels = fillVesselHoles(theseVessels, mask, tunnelCoord);
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

end

function vessel = fillVesselHoles(vessel, mask, tunnelCoord)%
% Dataset specific changes needed here, this function takes mask and
% drill holes into 'outer hull' to be able to use imfill, change drill
% location according to dataset :)
vesselMasked = or(vessel,mask);
if all(vesselMasked(1 : 80, tunnelCoord))
    warning('Tunnel doesn''t cut through');
end
vesselMasked(1 : 80, tunnelCoord) = 0;
vesselMasked(2080 : 2148, tunnelCoord) = 0;

vessel = imfill(vesselMasked, 'holes');% & ~mask;
end