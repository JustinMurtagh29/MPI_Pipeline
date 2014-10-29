function [parameter v] = setKDBparameter();

% Where to put data
parameter.writeLocation = '/nfs/bmo/mberning/knowledgeDB5/';
% Load global counter
parameter.counterLocation = '/zdata/manuel/data/cortex/kdbCounter.mat';

% Settings for dataset (written to knowledge DB)
parameter.settings.name = '2012-09-28_ex145_07x2';
parameter.settings.scale = [11.24 11.24 28];
parameter.settings.priority = 0;

% Border around problem to be looked at
% has to be big enough for rendering length calculation, sqrt(2)*maxLength in xy voxel equiv
parameter.border = [225 225 90];

% Viewports
v(1).name = 'tracer';
v(1).width = 200;
v(1).height = 200;
v(1).isRotated = 0;
v(2).name = 'b4b';
v(2).width = 312;
v(2).height = 214;
v(2).isRotated = 1;
v(2).openingAngle = 20;

end

