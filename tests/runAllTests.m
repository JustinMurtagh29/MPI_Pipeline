function runAllTests()
    % runAllTests
    %   Runs all tests in this directory.
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    % find current directory
    curDir = fileparts(mfilename('fullpath'));
    
    try
        % run all tests
        results = runtests(curDir, 'Recursively', true);
        display(results);
    catch e
        disp(getReport(e, 'extended'));
        exit(1);
    end
    
    % check if there are failed tasks
    isFailed = any(arrayfun( ...
        @(t) t.Failed > 0, results));
    exit(isFailed);
end