function param = getAvailibleParameterSettings(rootFolder)
% Pass cell array of folder to search for parameter settings

    % Initalize empty parameter struct
    nrParamsFound = 0;
    if nargin < 1
        rootFolder = {'/gaba/u/mberning/results/pipeline/' '/gaba/u/mberning/backup/zdata_manuel/results/pipeline/'};
    end
    for rF=1:length(rootFolder);
        files = dir(rootFolder{rF});
        for i=1:length(files)
            if exist([rootFolder{rF} files(i).name filesep 'allParameter.mat'], 'file')
                load([rootFolder{rF} files(i).name filesep 'allParameter.mat']);
                nrParamsFound = nrParamsFound + 1;
                param(nrParamsFound).p = changeFileLocationOfParameter(p, [rootFolder{rF} files(i).name filesep]);
                param(nrParamsFound).pT = changeFileLocationOfParameter(pT, [rootFolder{rF} files(i).name filesep]);
            end
        end
    end
    display(['Found ' num2str(nrParamsFound) ' pipeline parameter settings.']);
end

