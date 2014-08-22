function [p pT] = pipelineTest();
% Try implementing a run of the whole pipeline for some parameter settings here
% This is done in the logical order (could be optimized for time, e.g. start classification both in the beginning etc...)

% Create and set new parameters
[p, pT] = setParametersSettings();

% Pipeline on training regions 
j(1) = minicubeFwdPass(pT);
waitForState(j(1), 'finished');
j(2) = miniSegmentation(pT);
waitForState(j(2), 'finished');
j(3) = graphConstruction(pT);
waitForState(j(3), 'finished');
j(4) = miniFeature(pT);
waitForState(j(4), 'finished');
% Collect data of different training regions & normalize to [0 1]
prepareTrainingData(pT);
% Use data from graphs on training regions to initalize GP
prepareHyperparameter(pT);

% Same on test regions 
j(5) = minicubeFwdPass(p);
waitForState(j(5), 'finished');
j(6) = miniSegmentation(p);
waitForState(j(6), 'finished');
j(7) = graphConstruction(p);
waitForState(j(7), 'finished');
j(8) = miniFeature(p);
j(9) = correspondenceFinder(p);
waitForState(j(8), 'finished');
waitForState(j(9), 'finished');

% Start with some supervoxel graph stuff


end

