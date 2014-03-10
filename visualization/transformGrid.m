function [coords, coords3, coordsOld, matrix] = transformGrid( direction, position )
direction = direction .* [11.24 11.24 28];
direction = direction ./ norm(direction);
[X, Y, Z] = ndgrid(-2:2,-2:2,0:4);
coordsOld = [X(:) Y(:) Z(:)];

% ONB
orth1 = [1 0 -1].*direction([3 2 1]);
orth1 = orth1 ./ norm(orth1);
orth2 = cross(direction , orth1);
A = [orth1; orth2; direction];

% Rotation
coords = (A*coordsOld')';%coords/A;
% Skalierung
coords = coords .* repmat([1 1 11.24/28], size(coords, 1),1);
% Translation
coords = coords + repmat(position, size(coords,1),1);

% homogene Transformation
matrix = A;
matrix(3,:) = matrix(3,:) * 11.24/28;
matrix(1:3,4) = position;
matrix(4,1:3) = 0;
matrix(4,4) = 1;

coords2 = coordsOld;
coords2(:,4) = 1;
coords3 = (matrix* coords2')';

if all(all(abs(coords - coords3(:,1:3)) < 1e-3))
    display('success');
end

end