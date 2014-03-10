function rotatedVec = rotate(vector,axis,phi)
% vector 1 x 3, axis 1 x 3, not checked in here
L0=axis/norm(axis);
cphi=cosd(phi);
rotatedVec=vector*cphi+(vector*L0')*(1-cphi)*L0+cross(L0,vector)*sind(phi);
end

