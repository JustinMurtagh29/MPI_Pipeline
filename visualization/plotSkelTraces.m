function plotSkelTraces( in )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
colors = {'r' 'g' 'b' 'c' 'm' 'y' 'k'};
sizeLimit = size(in);
k = 1;
for i = 3:sizeLimit(1)
    for j = 3:sizeLimit(2)
        if ~isempty(in(i,j).trace)
            plotSkelTrace(in(i,j).trace{1}, colors{mod(k-1,7)+1});
            k = k + 1;
        end
    end
end
end

