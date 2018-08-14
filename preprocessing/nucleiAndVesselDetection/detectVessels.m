function vesselsPost = detectVessels( raw, visualize, tunnelCoord, threshold)
% Pass raw data array (tested with 07x2 mag4) and detect vessel in dataset,
% pass visualize flag in order to tune parameter on new dataset, goal
% should be to detect vessels with 0 FP rate, some FN are no problem, due to
% closing step in postprocessing

% By default: visualization off
if nargin < 2 || isempty(visualize)
    visualize = false;
end
if nargin < 3 || isempty(tunnelCoord)
    tunnelCoord = [];
end
if nargin < 4 || isempty(threshold)
    threshold = 150;
end

vessels = detectVesselsPre( raw, visualize, tunnelCoord, threshold);
vesselsPost = detectVesselsPost(vessels);