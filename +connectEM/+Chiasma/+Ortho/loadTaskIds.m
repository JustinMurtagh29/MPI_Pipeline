function tasks = loadTaskIds(taskFile)
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
end