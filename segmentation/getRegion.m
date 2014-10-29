function [x, y, z] = getRegion( centerOfMass2 )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
x=centerOfMass2(1)-6: centerOfMass2(1)+6;
y=centerOfMass2(2)-6: centerOfMass2(2)+6;
z=centerOfMass2(3)-3: centerOfMass2(3)+3;
x(x<1) = 1;
y(y<1) = 1;
z(z<1) = 1;
x(x>300) = 300;
y(y>300) = 300;
z(z>300) = 300;

end

