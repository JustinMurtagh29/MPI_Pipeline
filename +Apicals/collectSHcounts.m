% collectSHcounts 
% collect the spine heads counts for each agglomerate (designed for L4)
% -------------------------------------------------------------------------
% Inputs:
% agglos - agglomerates as sets of segments
% spineHeads - struct containing the spine head attachment by Alessandro
% Outputs:
% shCounts - Nx2 column1 agglo ind, column2 sh attached count
% -------------------------------------------------------------------------
% Author: Matej Zecevic <matej.zecevic@brain.mpg.de>
function [ shCounts ] = collectSHcounts( agglos, spineHeads )

% to get for all (small ones too) agglos = data.dendrites;
shCounts = zeros(length(agglos),1);
shCounts(1:end,1) = 1:length(agglos);
for d=1:length(agglos)
    shCounts(d) = numel(find(spineHeads.attached == d));
end

tic;
disp('collected spine head counts for agglos.'); toc;
disp('---------------------------------------');

end

