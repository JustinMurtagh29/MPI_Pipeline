function flights = filterFlights( flights )
%FILTERFLIGHTS Filter flights for patching into axons.
% INPUT flights: struct
%           Flights containing overlaps with agglos.
%           (see output of connectEM.Flight.overlapWithAgglos).
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>
% Based on code by Alessandro Motta

% remove flights with no or multiple attachments
numOverlaps = cellfun(@numel, flights.overlaps(:, 2));
Util.log('%5.1f %% of flights are dangling', 100 * mean(~numOverlaps));
Util.log('%5.1f %% of flights have too many overlaps', ...
    100 * mean(numOverlaps > 1));

flights = structfun( ...
    @(f) f(numOverlaps < 2, :), ...
    flights, 'UniformOutput', false);

% fill in zero for dangling flights
mask = cellfun(@isempty, flights.overlaps(:, 2));
flights.overlaps(mask, 2) = {0};
clear mask;

flights.overlaps = cell2mat(flights.overlaps);

% self-attachment is forbidden
mask = diff(flights.overlaps, 1, 2) ~= 0;
flights = structfun(@(f) f(mask, :), flights, 'UniformOutput', false);
clear mask;

% drop duplicate connections
[~, uniIds] = unique(sort(flights.overlaps, 2), 'rows', 'stable');
flights = structfun(@(f) f(uniIds, :), flights, 'UniformOutput', false);
clear uniIds;

end

