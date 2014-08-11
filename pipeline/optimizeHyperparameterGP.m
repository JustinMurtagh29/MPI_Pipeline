function optimizeHyperparameterGP(pT,mode)

if strcmp(mode,'edges')
    nrIndPoints = 1000;
    groundTrouthFile = pT.gp.initalGroundTruth;
    hypFile = pT.gp.hyperParameter;
elseif strcmp(mode,'glia')  
    nrIndPoints = 200;
    groundTrouthFile = pT.glia.initalGroundTruth;
    hypFile = pT.glia.hyperParameter;
end


% gpml toolbox (Rasmussen) startup script
run('/zdata/manuel/code/active/gpml/startup.m');
% Load normalized training data
load(groundTrouthFile);
% Define inducing points for FITC approximation
nrTrainingCases = size(trainingLabels,1);
rndIdx = randperm(nrTrainingCases);
indPoints = trainingData(rndIdx(1:nrIndPoints),:);
% Define all parameters of GP: mean, covariance, likelihood and inference functions
meanfunc = @meanConst;
hyp.mean = 0;
covfunc = {@covFITC, {@covSum, {@covLINard, @covNoise}}, indPoints};  % Try some 
hyp.cov = [zeros(size(trainingData,2),1); log(0.1)];
%covfunc = {@covFITC, {@covSEard}, indPoints};
%hyp.cov = log([ones(size(trainingData,2),1); 1]);
likfunc = @likErf;
inffunc = @infFITC_EP;

hyp = minimize(hyp, @gp, -100, inffunc, meanfunc, covfunc, likfunc, trainingData, trainingLabels);
save(hypFile, 'hyp', 'inffunc', 'meanfunc', 'covfunc', 'likfunc');

end

