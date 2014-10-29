function [x, y, z] = resizeBBox(x, y, z, sizeInput, minSize)

minSize = minSize - [length(x) length(y) length(z)];
minSize(minSize<0) = 0;
minSize = ceil(minSize/2);
x = [x(1)-minSize(1):1:x(1)-1 x x(end)+1:1:x(end)+minSize(1)];
y = [y(1)-minSize(2):1:y(1)-1 y y(end)+1:1:y(end)+minSize(2)];
z = [z(1)-minSize(3):1:z(1)-1 z z(end)+1:1:z(end)+minSize(3)];
x(x<1|x>sizeInput(1)) = [];
y(y<1|y>sizeInput(2)) = [];
z(z<1|z>sizeInput(3)) = [];

end

