function vesselsPost = detectVessels( raw, visualize, tunnelCoord, threshold)
% Pass raw data array (tested with 07x2 mag4) and detect vessel in dataset,
% pass visualize flag in order to tune parameter on new dataset, goal
% should be to detect vessels with 0 FP rate, some FN are no problem, due to
% closing step in postprocessing

% By default: visualization off
if nargin == 1
    visualize = false;
    tunnelCoord = 200;
end
if nargin <= 2
    tunnelCoord = 200;
end
if nargin <= 3
    threshold = 150;
end
vessels = detectVesselsPre( raw, visualize, tunnelCoord, threshold);
vesselsPost = detectVesselsPost(vessels);