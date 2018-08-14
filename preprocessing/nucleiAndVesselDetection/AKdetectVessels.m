function vesselsPost = AKdetectVessels( raw, visualize,output)
% Pass raw data array (tested with 07x2 mag4) and detect vessel in dataset,
% pass visualize flag in order to tune parameter on new dataset, goal
% should be to detect vessels with 0 FP rate, some FN are no problem, due to
% closing step in postprocessing

% By default: visualization off

vessels = AKdetectVesselsPre(raw, visualize);
save(fullfile(output,'vesselsPre'),'vessels','-v7.3');

vesselsPost = detectVesselsPost(vessels);