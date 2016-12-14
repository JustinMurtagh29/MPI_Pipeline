function spineHeadDetectionLocal(p,cubeNo)

% Please Note:
% This spine head classifier needs features calculated on segments / segment weights i.e. segmentWeightFile (instead of weightFile). You can get this segmentWeightFile
% by re-running calcFeatures and uncommenting segmentWeights parts 

me = mfilename;                                            % what is my filename
mydir = which(me); mydir = mydir(1:end-2-numel(me));        % where am I located

%Load classifier trained on 07x2 L4 dataset
load([mydir 'shdEnsemble.mat']);

% load the segment based features
m = load(p.local(cubeNo).segmentWeightFile);
segmentWeights = m.segmentWeights;

[~,scores]= predict(Ensemble,segmentWeights);

m=load(p.local(cubeNo).segmentFile);
segmentIds = [m.segments.Id]';

spineHeadScores = cat(2,double(segmentIds),scores);
Util.save([p.local(cubeNo).saveFolder 'spineHeadScores.mat'],spineHeadScores);

end
