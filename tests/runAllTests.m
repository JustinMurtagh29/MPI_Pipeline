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
    
    exit(0);
end