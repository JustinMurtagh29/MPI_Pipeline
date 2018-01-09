function [flights, wasDel] = rmEmptyAndCommentedFlights( flights )
%RMEMPTYANDCOMMENTEDFLIGHTS Remove empty and commented flights from flights
%struct.
% INPUT flights: struct
%           See output of connectEM.QueryUtil.loadProject
% OUTPUT flights: struct
%           The updated flights struct with flight paths removed that are
%           empty or contain a comment.
%       wasDel: [Nx1] logical
%           Logical array of the tasks that were deleted.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

wasDel = cellfun(@isempty, flights.nodes) | ...
     ~cellfun(@isempty, flights.comments);
flights = structfun(@(x)x(~wasDel), flights, 'uni', 0);

end

