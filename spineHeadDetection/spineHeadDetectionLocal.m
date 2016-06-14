function spineHeadDetectionLocal(p,cubeNo,Ensemble)

% Please Note:
% This spine head classifier needs features calculated on segments / segment weights i.e. segmentWeightFile (instead of weightFile). You can get this segmentWeightFile
% by re-running calcFeatures and uncommenting segmentWeights parts 

% load the wright file / features
m = load(p.local(cubeNo).segmentWeightFile);
segmentWeights = m.segmentWeights;

[~,scores]= predict(Ensemble,segmentWeights);

m=load(p.local(cubeNo).segmentFile);
segmentIds = [m.segments.Id]';

spineHeadScores = cat(2,double(segmentIds),scores);
save([p.local(cubeNo).saveFolder 'spineHeadScores.mat'],'spineHeadScores');

end
