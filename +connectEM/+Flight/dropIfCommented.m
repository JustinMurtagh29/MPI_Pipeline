function [flights, mask] = dropIfCommented(flights)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    mask = cellfun(@isempty, flights.comments);
    flights = structfun(@(f) f(mask), flights, 'UniformOutput', false);
end