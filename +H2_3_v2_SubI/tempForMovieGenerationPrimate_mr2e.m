
load('/tmpscratch/sahilloo/data/H2_3_v2_U1_SubI/pipelineRun_mr2e/allParameter.mat')
bbox = [7936, 8960, 1856,300,300,300];
%bbox = [8630, 7018, 1017, 200, 200, 200];
%outFiles = searchSegmentationThreshold(p,bbox,[0.7,0.11,0.18,0.26,0.30,0.35],false,true);
outFiles = searchSegmentationThreshold(p,bbox,[0.03,0.04],false,false);

