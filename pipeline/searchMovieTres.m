function searchMovieTres(Trange)
% function that allows you to create multiple movies of different
% segmentations. Set the appropriate Threshold in the Threshold loop. To
% prevent the multiple calculation of the classification a flag was
% introduced in makeSegmentationPreviewMovie


startup % starts everything and initiated Cluster Object

run configuration.m

flag = 0; % usually start with 1 whenever FOV was changed;
for threshold = Trange
    
    % change threshold and threshold function
    p.seg.threshold = threshold;
    p.seg.func = @(x)watershedSeg_v1_cortex(x,{p.seg.threshold 10});
    
    makeSegmentationPreviewMovie(p, [3933,4652; 994,2273 ; 1900 1999], flag) % second oint hast to be larger in every dimension
    
    flag = 0; % set flag to 0 to prevent recalculation of classification each time
end
end