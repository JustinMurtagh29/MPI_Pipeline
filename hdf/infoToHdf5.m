function infoToHdf5(file, info)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    isDirtyStrs = {'', ' (dirty)'};

    infoStr = cell(0, 1);
    infoStr{end + 1} = sprintf('%s', info.filename);

    for curRepoId = 1:numel(info.git_repos)
        curRepo = info.git_repos{curRepoId};
        infoStr{end + 1} = sprintf( ...
            '%s %s%s', curRepo.remote, curRepo.hash, ...
            isDirtyStrs{2 - isempty(curRepo.diff)}); %#ok
    end

    infoStr{end + 1} = sprintf( ...
        '%s@%s. MATLAB %s. %s', ...
        info.user, info.hostname, ...
        info.matlab_version, info.time);
    infoStr = strjoin(infoStr, newline);
    
    h5writeatt(file, '/', 'info', infoStr);
end
