function vessels = AKdetectVesselsPre(raw, visualize)
sizeRaw = size(raw);
vessels = false(size(raw));

if visualize
    figure('Position', [1921 1 1920 1124]);
    % Create all green RGB image for overlay
    green = cat(3, zeros([sizeRaw(1) sizeRaw(2)]), ones([sizeRaw(1) sizeRaw(2)]), zeros([sizeRaw(1) sizeRaw(2)]));
end

tic;
         
           threshold=160;
for z=1:sizeRaw(3)
    
%     if z<150
%         threshold=190;
%     else
%         if z<575
%         threshold=170;
%         else
%             if  z<630
%         threshold=160;
%             else
%         threshold=170;
%             end
%         end
%     end

    display(['Processing slice: ' num2str(z)]);
    thisImage = raw(:,:,z);
    theseVessels = false(size(thisImage));
    % Calculate mask (zeros voxel touching border of matrix)
    mask = ~imfill(~imdilate(thisImage == 0, ones([10 10])), 'holes');
    % Median Filter on image & keep only high pixel values
    thisImage = medfilt2(thisImage, [3 3]);
    thisImage(~mask & thisImage < threshold) = threshold;
    % These are only 129detection parameter that may need tuning on new dataset or
    % different magnification level (in addition to parameter in line above)
    regions = detectMSERFeatures(thisImage, 'ThresholdDelta', 1.5, 'MaxAreaVariation', 0.65, 'RegionAreaRange', [1000 300000]);
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

end

function vessel = fillVesselHoles(vessel, mask)
% Dataset specific changes needed here, this function takes mask and
% drill holes into 'outer hull' to be able to use imfill, change drill
% location according to dataset :)
vesselMasked = or(vessel,mask);
if all(vesselMasked(1 : 400, 560))
    warning('Tunnel doesn''t cut through');
end
vesselMasked(end-1400 : end, 560) = 0;
vesselMasked(1 : 80, 560) = 0;
vessel = imfill(vesselMasked, 'holes') & ~mask;
end