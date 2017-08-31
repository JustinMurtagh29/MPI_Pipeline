% filterSiteCriterion 
% filters all agglomerates that touch x boundaries (direction of pia and
% other side) and don't make contact to other boundaries (designed for L4)
% -------------------------------------------------------------------------
% Inputs:
% agglos - agglomerates as sets of segments
% points - agglomerates as points from segment meta
% bounds - 3x2 boundaries of dataset
% tol - tolerance around the border, the higher the more taken into account 
%       (in pixel so roughly 100 ~ 1 micron)
% Outputs:
% filtered agglos & points, rInd - indeces of filtered
% -------------------------------------------------------------------------
% Author: Matej Zecevic <matej.zecevic@brain.mpg.de>
function [ filteredAgglos, filteredPoints, rInd ] = filterSiteCriterion( agglos, points, bounds, tol )

tic;
% extract boundaries (contact sites)
xbounds = bounds(1,:); ybounds = bounds(2,:); zbounds = bounds(3,:);

% filter out agglos that touch the border in x direction
mask = zeros(length(agglos), 1);

for i=1:length(agglos)
    for j=1:size(points{i},1)
        point = points{i}(j);   % only consider x
        % if border is touched at both places within an agglo
        if point <= (xbounds(1) + tol)
            p1 = points{i}(j,:);
            for k=1:size(points{i},1)
                point = points{i}(k);
                p2 = points{i}(k,:);
                if point >= (xbounds(2) - tol)
                    if ~Apicals.moreContactPoints(points{i}, vertcat(p1,p2), tol, 0, tol, xbounds, ybounds, zbounds)
                        mask(i) = 1;
                    end
                end
            end
        end
    end
end
f1 = sum(mask) / length(agglos);
fprintf('filtered agglos, fraction of total: %f (%d/%d)\n', f1, sum(mask), length(agglos));

filteredAgglos = agglos(find(mask == 1));
filteredPoints = points(find(mask == 1));
rInd = find(mask ~= 0);

disp('finished execution of filterSiteCriterion.');
toc; disp('---------------------------------------');


end

