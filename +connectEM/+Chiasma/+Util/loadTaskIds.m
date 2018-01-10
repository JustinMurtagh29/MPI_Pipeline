function [tasks, startPos] = loadTaskIds(taskFile)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    fid = fopen(taskFile, 'rt');
    data = fread(fid, 'char=>char');
    fclose(fid);
    
    % split into lines
    data = reshape(data, 1, []);
    data = strsplit(data, '\n');
    
    % remove header
    data(1) = [];
    
    % split into fields
    data = cellfun( ...
        @(s) strsplit(s, ','), ...
        data, 'UniformOutput', false);
    data = cat(1, data{:});
    
    % extract first column with task IDs
    tasks = table;
    tasks.id = data(:, 1);
    tasks.nmlFile = data(:, 2);
    
    if nargout > 1 % read out start position
        
        % get rid of commas
        data(:,3) = cellfun(@(x)x(2:end), data(:,3), 'uni', 0);
        data(:,5) = cellfun(@(x)x(1:strfind(x, ')')-1), data(:,5), 'uni', 0);
        
        % to array
        startPos = cell2mat(cellfun(@str2num, data(:,3:5), 'uni', 0));
    end
end