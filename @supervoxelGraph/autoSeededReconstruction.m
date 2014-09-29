function seededReconstruction(p, adjMatrix, startSv)
%not finished yet
th = 0.4;    
reconstruction = struct('propability', [], 'components', []);
index = 1;
[p, idx] = max(adjMatrix(startSv, :));
prob = [];
x = startSv;
y = idx;

while adjMatrix(x,y) >= th
    reconstruction(index).propability = adjMatrix(x,y);
    prob = [prob; adjMatrix(x,y)];
    reconstruction(index).components = [x,y];

    if mod(index,2) == 1
        [p, idx] = max(adjMatrix([1:x-1 x+1:end], y));
        x = idx;
        y = y;
    else
        [p, idx] = max(adjMatrix(y,[1:x-1 x+1:end]));
        x = y;
        y = idx;
    end
    index = index + 1;
end

end
