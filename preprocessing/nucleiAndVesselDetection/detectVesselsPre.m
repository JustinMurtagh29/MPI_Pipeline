function vessels = detectVesselsPre(raw, visualize, tunnelCoord, threshold,mag)
sizeRaw = size(raw);
vessels = false(size(raw));
if ~exist('mag','var') || isempty(mag)
    mag = 1;
end
if visualize
    figure;%('Position', [3841 1 1920 999]);
    % Create all green RGB image for overlay
end
green = cat(3, zeros([sizeRaw(1) sizeRaw(2)]), ones([sizeRaw(1) sizeRaw(2)]), zeros([sizeRaw(1) sizeRaw(2)]));
tic;
parfor z=1:sizeRaw(3)
    display(['Processing slice: ' num2str(z)]);
    thisImage = raw(:,:,z);
    theseVessels = false(size(thisImage));
    % Calculate mask (zeros voxel touching border of matrix) -> one could
    % improve this if the mask is generated during cubing, as it would
    % anyways be good to always know what is outside the dataset...
    mask = ~imfill(~imdilate(thisImage == 0, ones([10 10])), 'holes'); 
    % Median Filter on image & keep only high pixel values
    thisImage = medfilt2(thisImage, [3 3]);
    
    % These are only detection parameter that may need tuning on new dataset or
    % different magnification level (in addition to parameter in line above)
    regions = detectMSERFeatures(thisImage, 'ThresholdDelta', 2, 'MaxAreaVariation', 0.5, 'RegionAreaRange', [1000*16 1e6*16]/mag^2);
    
%     figure;histogram(cellfun(@(x) mean(thisImage(sub2ind(size(thisImage), x(:,2), x(:,1)))),regions.PixelList),0:5:255);
    indHigh = cellfun(@(x) mean(thisImage(sub2ind(size(thisImage), x(:,2), x(:,1)))),regions.PixelList) > threshold;
    regions = regions(indHigh);
    % Set regions to true that were detected by MSER
    for i=1:regions.Count
        idx = sub2ind(size(thisImage), regions.PixelList{i}(:,2), regions.PixelList{i}(:,1));
        theseVessels(idx) = 1;
    end
    % Fill holes in detection
    theseVessels = fillHoles(theseVessels, mask, tunnelCoord);
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