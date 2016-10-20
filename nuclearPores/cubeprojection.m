function [Nver, Nhor] = cubeprojection(nuclei_roi, raw_roi)
%Triangulate isosurface
p = isosurface(nuclei_roi,.5);
TR=triangulation(p.faces,p.vertices);
U = incenter(TR);
fn = faceNormal(TR);
%%
delta = 4;
a = zeros(size(U,1),delta);
for i = 1:delta
    u = floor(U + (i-2)*fn);
    a(:,i) = raw_roi(sub2ind(size(raw_roi),u(:,2),u(:,1),u(:,3)));
end
a = mean(a,2);

extend = range(U)/max(range(U));
center = mean(U);
P = [(U(:,1)-center(1))/extend(1),(U(:,2)-center(2))/extend(2),(U(:,3)-center(3))/extend(3)];
[thetaVer,phiVer] = cart2sph(-P(:,3),P(:,1),P(:,2)); %polar on X-Z plane, rotate nucleus along Y axis
[thetaHor,phiHor] = cart2sph(P(:,2),P(:,1),P(:,3)); %polar on X-Y plane
NprojectV = getprojection(thetaVer, phiVer, a);
NprojectH = getprojection(thetaHor, phiHor, a);
Nver = [NprojectV(:,3151:3600),NprojectV(:,1:2250)];
Nhor = [NprojectH(:,1351:3600),NprojectH(:,1:450)];
