function devStepsPipeline( dataFolder, saveFolder, steps )
%DEVSTEPSPIPELINE Train all classifiers for dev steps.
% INPUT dataFolder: string
%           Path to training data folder.
%       saveFolder: string
%           Path to folder where all data will be saved.
%       steps: (Optional) [1xN] int
%           Specify which dev steps to calculate.
%           (Default: 1:8)
%
% NOTE The feature map must be defined in Paper.SynEM.devStepsFM.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

saveFolder = SynEM.Util.addFilesep(saveFolder);

if ~exist('steps', 'var') || isempty(steps)
    steps = 1:8;
end

if iscolumn(steps); steps = steps'; end

if ~exist('saveFolder','dir')
    mkdir(saveFolder);
end

for s = steps
    fprintf('[%s] Processing step %d.\n', datestr(now), s);
    switch s
        case {1, 2, 3, 4, 5}
            labelType = 'single';
            ensembleArgs = {'method','AdaBoostM1'};
        case 6
            labelType = 'direction';
            ensembleArgs = {'method','LogitBoost'};
        case 7
            labelType = 'augment';
            ensembleArgs = {'method','LogitBoost'};
        case 8
            labelType = 'single';
            ensembleArgs = {'method','LogitBoost'};
        otherwise
            labelType = [];
            ensembleArgs = [];
    end
    fm = Paper.SynEM.devStepsFM(s);
    %label mode is currently selected when loading the label (and always
    %requires 'direction'), i.e. only depending on labelType
    fm.mode = 'direction';
    currSaveFolder = [saveFolder 'DevStep' num2str(s)];
    SynEM.Paper.trainingPipeline(dataFolder, currSaveFolder, fm, true, ...
        labelType, ensembleArgs);
end



end

