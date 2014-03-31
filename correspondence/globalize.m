function result = globalize(cubeCoords1, cubeCoords2, cubeSize, overlap, seg_path )
%GLOBALIZE 
% Create index array of overlapping segments

cube1 = load([seg_path,'x000',num2str(cubeCoords1(1)),'y000',num2str(cubeCoords1(2)),'z000',num2str(cubeCoords1(3)),'\seg.mat']);
cube1 = cube1.seg;
cube2 = load([seg_path,'x000',num2str(cubeCoords2(1)),'y000',num2str(cubeCoords2(2)),'z000',num2str(cubeCoords2(3)),'\seg.mat']);
cube2 = cube2.seg;
coordShift = cubeCoords2 - cubeCoords1;

% Calculating coordinates for overlap1 and overlap 2
for i = 1:3 
    x1(i,1) = cubeSize(i)/4 + 0.5*cubeSize(i)*any(coordShift(i));
    x1(i,2) = cubeSize(i)*3/4 + any(coordShift(i));
    
    x2(i,1) = cubeSize(i)/4 - any(coordShift(i));
    x2(i,2) = cubeSize(i)*3/4 - any(coordShift(i))*0.5*cubeSize(i);
end

overlap1 = cube1(x1(1,1) : x1(1,2), x1(2,1) : x1(2,2), x1(3,1) : x1(3,2));
overlap2 = cube2(x2(1,1) : x2(1,2), x2(2,1) : x2(2,2), x2(3,1) : x2(3,2));

%create one dimensional array to compare both overlaps
[x,y,z] = size(overlap1);
ind1 = [];
ind2 = [];
for x = 1:x
    for y = 1:y
        for z = 1:z
            ind1 = [ind1; overlap1(x,y,z)];
            ind2 = [ind2; overlap2(x,y,z)];
        end
    end
end

%create output struct
ind = [ind1 ind2];
[ind(:,1), idx] = sort(ind(:,1));
ind(:,2) = ind(idx,2);
ind(any(ind == 0,2),:) = [];
result.ind = uint32(unique(ind, 'rows'));
result.ind(:,1) = result.ind(:,1) + 10000000 * cubeCoords1(1) +  1000000 * cubeCoords1(2) +  100000 * cubeCoords1(3); 
result.ind(:,2) = result.ind(:,2) + 10000000 * cubeCoords2(1) +  1000000 * cubeCoords2(2) +  100000 * cubeCoords2(3);

result.cubeCoords1 = cubeCoords1;
result.cubeCoords2 = cubeCoords2;

result.overlap1 = overlap1;
result.overlap2 = overlap2;
