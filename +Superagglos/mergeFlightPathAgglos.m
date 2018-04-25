function axons = mergeFlightPathAgglos( ov, axonsBig, axonsSmall, toAgglo )
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
%       toAgglo: (Optional) logical
%           Flag indicating that the output should be in the agglo format,
%           i.e. a list of segmentation ids for each agglo. This is
%           typically much faster.
%           Note that Superagglos.mergeOnOverlaps results in different
%           output agglos (should currently not be used anymore).
%           (Default: true)
% OUTPUT axons: [Nx1] struct
%           Resulting axons in the superagglo format.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if ~exist('toAgglo', 'var') || isempty(toAgglo)
    toAgglo = true;
end

axons = axonsBig;
for i = 1:length(axonsBig)
    if toAgglo
        if ~isempty(ov{i}) % simply concatenate nodes
            % assumes axons seg ids are mutually exclusive
            axons(i).nodes = cat(1, axons(i).nodes, ...
                cell2mat({axonsSmall(ov{i}).nodes}'));
        end
    else
        for j = 1:length(ov{i})
            % axonsBig are used as second input such that the flight path
            % from the big axon is removed and replaced with the small axon
            axonsBig(i) = Superagglos.mergeOnOverlaps( ...
                axonsSmall(ov{i}(j)), axonsBig(i), ...
                'scale', [11.24, 11.24, 28], ...
                'overlapDistNm', 500, ...
                'use2015b', true);
        end
    end
end

if toAgglo
    axons = Superagglos.getSegIds(axons);
end

end

