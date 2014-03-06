function [ result ] = getOverlaps( cubeLims )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

overlaps = [];
sizeDims = [cubeLims(1,2)-cubeLims(1,1)+1 cubeLims(2,2)-cubeLims(2,1)+1 cubeLims(3,2)-cubeLims(3,1)+1];

%i did sub2ind and later ind2sub so i could do unique on the rows (i didnt
%know how to with 6 values..:P)
for x = cubeLims(1,1):cubeLims(1,2)
    for y = cubeLims(2,1):cubeLims(2,2)
        for z = cubeLims(3,1):cubeLims(3,2)
            coords = [x y z];
            ind = sub2ind(sizeDims,x,y,z);
            if(coords(1) - 1 >= cubeLims(1,1))
                overlaps = [overlaps; ind sub2ind(sizeDims,x-1,y,z)];
            end
            if(coords(1) + 1 <= cubeLims(1,2))
                overlaps = [overlaps; ind sub2ind(sizeDims,x+1,y,z)];
            end
            
            if(coords(2) - 1 >= cubeLims(2,1))
                overlaps = [overlaps; ind sub2ind(sizeDims,x,y-1,z)];
            end            
            if(coords(2) + 1 <= cubeLims(2,2))
                overlaps = [overlaps; ind sub2ind(sizeDims,x,y+1,z)];
            end
            
            if(coords(3) - 1 >= cubeLims(3,1))
                overlaps = [overlaps; ind sub2ind(sizeDims,x,y,z-1)];
            end
            if(coords(3) + 1 <= cubeLims(3,2))
                overlaps = [overlaps; ind sub2ind(sizeDims,x,y,z+1)];
            end
            
        end
    end
end

overlaps = unique(sort(overlaps,2),'rows');

result = [];
for i = 1 : size(overlaps,1)
    [x,y,z] = ind2sub(sizeDims, overlaps(i,1));
    [x2,y2,z2] = ind2sub(sizeDims, overlaps(i,2));
    result = [result; x y z x2 y2 z2]; 
end
        
