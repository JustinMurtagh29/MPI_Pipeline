function [axons, parentIds, ov] = ...
        mergeFlightPathAgglos( ov, axonsBig, axonsSmall, mergeMode )
%MERGEFLIGHTPATHAGGLOS Merge the small axon agglomerates that are picked up
% by a flight path into the large axons class.
% INPUT ov: [Nx1] cell
%           Linear indices of axonsSmall that overlap with the
%           corresponding axonBig.
%           see output of Superagglos.flightPathAggloOverlap
%       axonsBig: struct
%           The large axons in the superagglo format.
%       axonsSmall: struct
%           The small axons in the superagglo format that are merged into
%           the large axons.
%       mergeMode: (Optional) string
%           Specification of the merge mode:
%           'toAgglo': Only merges the segment ids which is typically
%               rather fast (Default)
%           'overlaps': see Superagglos.mergeOnOverlaps
%               Note that Superagglos.mergeOnOverlaps results in different
%               output agglos in terms of segment ids (should currently not
%               be used anymore).
%           'simple': Merges two superagglos without removing any 
%                redundancies. This can result in duplicate flight paths
%                etc.
% OUTPUT axons: [Nx1] struct
%           Resulting axons in the superagglo format.
%        parentIds: [Nx1] double
%           Linear indices of the `axonBig` corresponding to `axons`.
%        ov: [Nx1] cell
%           The updated overlap struct (in case flight path only agglos
%           were deleted when using toAgglo = true).
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if ~exist('mergeMode', 'var') || isempty(mergeMode)
    mergeMode = 'toAgglo';
end

axons = axonsBig;
parentIds = reshape(1:numel(axons), size(axons));


switch mergeMode
    
    case 'toAgglo'
        % simply concatenate nodes
        for i = 1:length(axonsBig)
            if ~isempty(ov{i})
                
                % assumes axons seg ids are mutually exclusive
                axons(i).nodes = cat(1, axons(i).nodes, ...
                    cell2mat({axonsSmall(ov{i}).nodes}'));
            end
        end
        toDel = Superagglos.getPureFlightPathAgglos(axons);
        axons(toDel) = [];
        parentIds(toDel) = [];
        ov(toDel) = [];
        axons = Superagglos.getSegIds(axons);
        
    case 'simple'
        % simply concatenate nodes & add random edge
        for i = 1:length(axonsBig)
            if ~isempty(ov{i})
                nodes = axons(i).nodes;
                % assumes axons seg ids are mutually exclusive
                axons(i).nodes = cat(1, axons(i).nodes, ...
                    cell2mat({axonsSmall(ov{i}).nodes}'));
                f = @(x) x(:, 1:3) .* [11.24, 11.24, 28];
                for j = 1:length(ov{i})
                    d = pdist2(f(nodes), f(axonsSmall(ov{i}(j)).nodes));
                    [~, idx] = min(d(:));
                    [newEdge(1), newEdge(2)] = ind2sub(size(d), idx);
                    offset = max(axons(i).edges(:));
                    newEdge(2) = newEdge(2) + offset;
                    axons(i).edges = cat(1, axons(i).edges, ...
                        axonsSmall(ov{i}(j)).edges + offset, ...
                        newEdge(:)');
                end
            end
        end
        
    case 'overlap'
        for i = 1:length(axonsBig)
            for j = 1:length(ov{i})
                % axonsBig are used as second input such that the flight path
                % from the big axon is removed and replaced with the small axon
                axonsBig(i) = Superagglos.mergeOnOverlaps( ...
                    axonsSmall(ov{i}(j)), axonsBig(i), ...
                    'scale', [11.24, 11.24, 28], ...
                    'overlapDistNm', 500);
            end
        end
    otherwise
        error('Unknown merge mode %s.', mergeMode);
end

end

