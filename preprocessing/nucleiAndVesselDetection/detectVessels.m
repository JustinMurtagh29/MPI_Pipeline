function vesselsPost = detectVessels( raw, threshold,mag,visualize, tunnelCoord)
% Pass raw data array (tested with 07x2 mag4) and detect vessel in dataset,
% pass visualize flag in order to tune parameter on new dataset, goal
% should be to detect vessels with 0 FP rate, some FN are no problem, due to
% closing step in postprocessing
%
% INPUT
% raw           raw data
% threshold     threshold above which mean region intensity has to be to be
%               detected as BV
% mag           magnification of raw, necessary for size-dependent steps
% visualize     by default: visualization off
% tunnelCoord   By default in the middle of the plane a tunnel is drilled
%               through the black borders in order to fill holes. If this
%               fails whyever, it can be given manually
%
% modified by marcel.beining@brain.mpg.de

if nargin < 2 || isempty(threshold)
    threshold = 220;
end
if nargin < 3 || isempty(mag)
    mag = 1;
end
if nargin < 4 || isempty(visualize)
    visualize = false;
end
if nargin < 5 || isempty(tunnelCoord)
    tunnelCoord = [];
end

vessels = detectVesselsPre( raw, visualize, tunnelCoord, threshold,mag);
vesselsPost = detectVesselsPost(vessels);