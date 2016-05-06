function [ classificationStruct ] = getPredictionVariables( classifier, method )
%GETPREDICTIONVARIABLES Extract the variables important for prediction.
% Extract the variables used for classification from an ensemble of
% decision stumps trained with fitensemble (currently for LogitBoost and
% AdaBoostM1) using decision stumps.
% Use the classificationStruct with the predict function to classify new
% datapoints.
% INPUT classifier: The classifier.
%       method: String which methods was used (currently only LogitBoost
%               and AdaBoost.

switch method
    case 'LogitBoost'
        T = classifier.NumTrained;
        W = classifier.TrainedWeights;
        cutVars = zeros(T,1);
        cutPoints = zeros(T,1);
        leftNodePred = zeros(T,1);
        rightNodePred = zeros(T,1);
        for t = 1:T
            weakLearner = classifier.Trained{t}.CompactRegressionLearner;
            cutPoints(t) = weakLearner.CutPoint(1);
            cutVars(t) = str2double(weakLearner.CutVar{1}(2:end));
            leftNodePred(t) = -weakLearner.NodeMean(2);
            rightNodePred(t) = -weakLearner.NodeMean(3);
        end
    case 'AdaBoostM1'
        T = classifier.NumTrained;
        W = classifier.TrainedWeights;
        cutVars = zeros(T,1);
        cutPoints = zeros(T,1);
        leftNodePred = zeros(T,1);
        rightNodePred = zeros(T,1);
        for t = 1:T
            weakLearner = classifier.Trained{t};
            cutPoints(t) = weakLearner.CutPoint(1);
            cutVars(t) = str2double(weakLearner.CutVar{1}(2:end));
            if strcmp(weakLearner.NodeClass{2},'true')
                leftNodePred(t) = 1;
            else
                leftNodePred(t) = -1;
            end
            if strcmp(weakLearner.NodeClass{3},'true')
                rightNodePred(t) = 1;
            else
                rightNodePred(t) = -1;
            end
        end
end

classificationStruct.W = W;
classificationStruct.cutVars = cutVars;
classificationStruct.cutPoints = cutPoints;
classificationStruct.leftNodePred = leftNodePred;
classificationStruct.rightNodePred = rightNodePred;

end

