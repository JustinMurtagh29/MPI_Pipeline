% script to generate the axon agglo version
% 'axons_19_a_partiallySplit_v3'
% which is simply a copy of axons_19_a_partiallySplit_v2 to redo the
% flight path pick-up
%
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>


info = Util.runInfo();

rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
axonFile = fullfile(rootDir, 'aggloState', ...
    'axons_19_a_partiallySplit_v2.mat');
outfile = fullfile(rootDir, 'aggloState', ...
    'axons_19_a_partiallySplit_v3.mat');
m = load(axonFile);
m.info = info;

save(outfile, '-struct', 'm');
Util.protect(outfile);

