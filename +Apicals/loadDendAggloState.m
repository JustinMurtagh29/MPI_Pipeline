% loadDendAggloState 
% loads all necessary states to run the filter (designed for L4)
% -------------------------------------------------------------------------
% Inputs:
% dendAggloStatePath - current state of dendrite agglomerates
% spineHeadsStatePath (optional) - current state of spine heads attachment
% dendLensPath - lengths for current state of dendrite agglos
% parameterPath - path to parameters file of dataset
% metaPath - path to segmentMeta file
% local - when using from mount point e.g. '/mnt/gaba'
% Outputs:
% all necessary data for the filtering 
% -------------------------------------------------------------------------
% Author: Matej Zecevic <matej.zecevic@brain.mpg.de>
function [ agglosNew, agglos, agglos_all, points, indBigDends, bounds, spineHeads, lens, sM, p ] = ...
    loadDendAggloState( dendAggloStatePath, spineHeadsStatePath, dendLensPath, parameterPath, metaPath, local )
tic;
% collect the dendrite agglos (file from Marcel in new format) only > 5 um
data = load(dendAggloStatePath);
%agglosNew_all = data.newAgglos; % for my file with removed soma
agglosNew_all = data.dendrites;
if exist('indBigDends','var')
    agglosNew = data.dendrites(data.indBigDends); 
    indBigDends = data.indBigDends;

<<<<<<< HEAD
    % generate all agglo representation for following code
    agglos = connectEM.transformAggloNewOldRepr(agglosNew);
    agglos_all = connectEM.transformAggloNewOldRepr(agglosNew_all);
else
    indBigDends = NaN;
    agglosNew = agglosNew_all;
    agglos = connectEM.transformAggloNewOldRepr(agglosNew_all);
    agglos_all = agglos;
end
=======
% generate all agglo representation for following code
agglos = Superagglos.transformAggloNewOldRepr(agglosNew);
agglos_all = Superagglos.transformAggloNewOldRepr(agglosNew_all);
>>>>>>> cd5baf8448c1ab8eb91287ace82a1fd3b3399f0e

% all lengths (number of seg ids) of the agglos
l = cellfun(@length, agglos);

if exist('spineHeadsStatePath','var')
    spineHeads = load(spineHeadsStatePath);
end

% collect parameter and bbox
if exist('parameterPath','var')
    p = load(parameterPath);
    p = p.p;
    if exist('local','var')
        p.saveFolder = [local p.saveFolder];
        p.seg.root = [local p.seg.root]; %'/gaba/wKcubes.new/Connectomics department/2012-09-28_ex145_07x2_ROI2017/segmentation/1/'
    end
    bounds = p.bbox;
end

% collect points from segments
if exist('metaPath','var')
    sM = load(metaPath);
    points = cellfun(@(x) sM.point(:,x)',agglos,'uni',0);
end

% collect dendrite agglo lengths (in um)
if exist('dendLensPath','var')
    load(dendLensPath);
    lens = lens / 1E3;
end

disp('finished execution of loadDendAggloState.');
toc; disp('---------------------------------------');

    
end

