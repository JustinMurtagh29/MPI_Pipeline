function contOptimizeHyperparameterGP(pT,mode)

global GLOBAL_CODE_DIR;

timerVal = tic;

if strcmp(mode,'edges')
  groundTrouthFile = pT.gp.initalGroundTruth;
  hypFile = pT.gp.hyperParameter;
elseif strcmp(mode,'glia')  
  groundTrouthFile = pT.glia.initalGroundTruth;
  hypFile = pT.glia.hyperParameter;
end


% gpml toolbox (Rasmussen) startup script
run([GLOBAL_CODE_DIR 'active/gpml/startup.m']);
% Load normalized training data
load(groundTrouthFile);
load(hypFile);

hyp = minimize(hyp, @gp, -20, inffunc, meanfunc, covfunc, likfunc, trainingData, trainingLabels);
save([hypFile(1:end-4) '-' datestr(clock, 30) '.mat'], 'hyp', 'inffunc', 'meanfunc', 'covfunc', 'likfunc');

toc(timerVal)

end

