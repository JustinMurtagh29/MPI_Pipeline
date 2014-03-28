function [ result ] = getOverlaps( cubeLims )
%UNTITLED 
%  find overlaps between neighbouring cubes

result = [];
for x = cubeLims(1,1):cubeLims(1,2)
    for y = cubeLims(2,1):cubeLims(2,2)
        for z = cubeLims(3,1):cubeLims(3,2)

            coords = [x y z];
            if(coords(1) + 1 <= cubeLims(1,2))
                result = [result; coords x+1 y z];
            end
                     
            if(coords(2) + 1 <= cubeLims(2,2))
                result = [result; coords x y+1 z];
            end
            
            if(coords(3) + 1 <= cubeLims(3,2))
               result = [result; coords x y z+1];
            end
            
        end
    end
end
