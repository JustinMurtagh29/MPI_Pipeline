function AlgComp = collectDevStepsResults( saveFolder, steps )
%COLLECTDEVSTEPSRESULTS Collect the results from the devSteps pipeline.
% INPUT saveFolder: string
%           See Paper.SynEM.devStepsPipeline.
%       steps: (Optional) [1xN] int
%           Specify which dev steps to calculate.
%           (Default: 1:8)
% OUTPUT AlgComp: [1xN] struct
%           Struct containing the fields
%           'classifier': The compact classifier.
%           'scores': Validation scores
%           'rp': Validation rp.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if ~exist('steps', 'var') || isempty(steps)
    steps = 1:8;
end
if iscolumn(steps); steps = steps'; end
saveFolder = SynEM.Util.addFilesep(saveFolder);

[~, yVal] = SynEM.Util.getTrainingDataFrom([saveFolder, 'DevStep' ...
    num2str(steps(1)) filesep 'val'], 'single', true);

AlgComp = struct;
for s = steps
    m = load([saveFolder, 'DevStep' num2str(s) filesep 'classifier' ...
        filesep 'LogitBoost.mat']);
    AlgComp(s).classifier = m.classifier;
    m = load([saveFolder, 'DevStep' num2str(s) filesep 'classifier' ...
        filesep 'result.mat']);
    AlgComp(s).scores = m.scoresVal;
    switch s
        case {1, 2, 3, 4, 5, 7 8}
            AlgComp(s).rp = SynEM.Eval.interfaceRP(cell2mat(yVal), ...
                AlgComp(s).scores);
        case 6 %direction
            AlgComp(s).scores = m.scoresVal;
            AlgComp(s).rp = SynEM.Eval.interfaceRPDir2Single(yVal, ...
                AlgComp(s).scores);
    end
end


end

