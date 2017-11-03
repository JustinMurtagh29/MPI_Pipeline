function outSkels = SkelOverlap(inputDir,superagglos,p)
% This function searches for each superagglo that overlaps the most with
% each of the skeletons (nml files, 1 skel per file) found in inputDir and
% writes them in the same folder as skeletons.
%
% INPUT
% inputDir      string with path to nml files which are used
% superagglos   agglos in the superagglo format which are searched through
% p             (optional) parameter struct of the pipeline run. If not
%               defined, the L4 20170217_ROI pipeline run is used
%
% OUTPUT
% outSkels  the superagglos which overlapped as skeletons
%
% By marcel.beining@brain.mpg.de


if ~exist('p','var') || isempty(p)
    load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat');
end
[superaggloLUT,superaggloSegIds] = Superagglos.buildLUT(superagglos);

if ~exist(inputDir,'dir')
    error('Folder %s is not existent',inputDir);
end
files = dir(fullfile(inputDir,'*.nml'));

outSkels = skeleton();
for f = 1:numel(files)
    skel = skeleton(fullfile(inputDir,files(f).name));  % load skeleton
    skelCoords = cell2mat(cellfun(@(x) x(:,1:3),skel.nodes,'uni',0));  % putting all skel nodes together
    warning('OFF','auxiliaryMethods:readKnossosCube')
    skelSegIds = unique(Seg.Global.getSegIds(p,skelCoords));  % extract the seg Ids of all skel nodes but use each segId only once
    warning('ON','auxiliaryMethods:readKnossosCube')
    
    [~,aggloOverlapsSkel] = ismember(skelSegIds,superaggloSegIds,'rows'); % check which superagglos (segIds) have overlaps with skeleton (segIds)
    ind = mode(superaggloLUT(superaggloSegIds(nonzeros(aggloOverlapsSkel)))); % get superagglo that overlaps the most
    if ind == 0
        error('No superagglo was found that has segID overlaps with the skeleton from file %s!',files(f).name)
    end
    
    % transform found superagglo to skeleton and write it to file
    skelSuperagglo = Superagglos.toSkel(superagglos(ind));
    outSkels = outSkels.addTreeFromSkel(skelSuperagglo);
    skelSuperagglo.write(fullfile(inputDir,strrep(files(f).name,'.nml',sprintf('_Superagglo_%d.nml',ind)) ));
end