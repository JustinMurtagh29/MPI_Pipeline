function [agglo,aggloLUT] = remove(agglo,ind,aggloLUT)
% remove an agglo from the list of agglos and its LUT


if nargin > 2 && nargout > 1
    aggloLUT = connectEM.changem(aggloLUT,((0:numel(agglo))-[0, cumsum(accumarray(ind+1,1,[numel(agglo),1]))']).*~[0;accumarray(ind,1,[numel(agglo),1])]',0:numel(agglo));
end
agglo(ind) = [];