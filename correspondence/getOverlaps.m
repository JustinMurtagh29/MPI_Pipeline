function result = getOverlaps( numberCubes )
% getOverlaps: Construct array indicating touching faces between segmentation cubes 

result = [];
for x = 1:numberCubes(1)
    for y = 1:numberCubes(2)
        for z = 1:numberCubes(3)
            coords = [x y z];
            if(coords(1) + 1 <= numberCubes(1))
                result = [result; coords x+1 y z];
            end
            if(coords(2) + 1 <= numberCubes(2))
                result = [result; coords x y+1 z];
            end
            if(coords(3) + 1 <= numberCubes(3))
               result = [result; coords x y z+1];
            end
        end
    end
end
