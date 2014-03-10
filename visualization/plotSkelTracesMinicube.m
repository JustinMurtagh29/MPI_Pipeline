function plotSkelTracesMinicube( in, posM, sizeM, color )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
sizeLimit = size(in);
for i = 1:sizeLimit(1)
    for j = 1:sizeLimit(2)
        if ~isempty(in(i,j).trace) 
            plotSkelTraceMinicube(in(i,j).trace{1}, posM, sizeM, color);
        end
    end
end
end

