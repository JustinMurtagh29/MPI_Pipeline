function setup()
%SETUP Setup operations for SynEM.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

% go to directory of this file
prevDir = pwd();
thisDir = fileparts(mfilename('fullpath'));
cd(thisDir);

% compile mex files
mex CXXFLAGS='$CXXFLAGS -std=c++11' ...
    -outdir +Aux...
    +Aux/eig3S.cpp
mex CXXFLAGS='$CXXFLAGS -std=c++11'...
    -outdir +Aux...
    +Aux/sortAbs.cpp

% go to original directory
cd(prevDir);


end

