function matrix4 = angle_fudge2(direction_vector)
% generate rotation matrix from direction vector
% discontinuous in [0, 0, z]
% Kevin M. Boergens, 2016

assert(norm(direction_vector) > 0)
matrix4 = [-1, 0, 0, 0, 0, 0, 0.4014285632542201, 0, 0, 1, 0, 0, 0, 0, 0, 1];
firstComp = -asin(direction_vector(3) / norm(direction_vector));
if norm(direction_vector(1 : 2))
    secondComp = feval(@(x)x(4)*x(3), vrrotvec([0, 1, 0], [direction_vector(1 : 2), 0]));
else
    secondComp = 0;
end

pitchVec = [1, 0, 0];
yawVec = [0, 1, 0];
matrix4 = M4x4_rotate(secondComp, yawVec, matrix4, matrix4);
matrix4 = M4x4_rotate(firstComp, pitchVec, matrix4, matrix4);
end

%adatped from webgl-mjs
function r = M4x4_rotate (angle, axis, m, r)
a0 = axis(1);
a1 = axis(2);
a2 = axis(3);
l = sqrt(a0 * a0 + a1 * a1 + a2 * a2);
x = a0;
y = a1;
z = a2;
if (l ~= 1.0)
    im = 1.0 / l;
    x = x * im;
    y = y * im;
    z = z * im;
end
c = cos(angle);
c1 = 1 - c;
s = sin(angle);
xs = x * s;
ys = y * s;
zs = z * s;
xyc1 = x * y * c1;
xzc1 = x * z * c1;
yzc1 = y * z * c1;

m11 = m(1);
m21 = m(2);
m31 = m(3);
m41 = m(4);
m12 = m(5);
m22 = m(6);
m32 = m(7);
m42 = m(8);
m13 = m(9);
m23 = m(10);
m33 = m(11);
m43 = m(12);

t11 = x * x * c1 + c;
t21 = xyc1 + zs;
t31 = xzc1 - ys;
t12 = xyc1 - zs;
t22 = y * y * c1 + c;
t32 = yzc1 + xs;
t13 = xzc1 + ys;
t23 = yzc1 - xs;
t33 = z * z * c1 + c;

r(1) = m11 * t11 + m12 * t21 + m13 * t31;
r(2) = m21 * t11 + m22 * t21 + m23 * t31;
r(3) = m31 * t11 + m32 * t21 + m33 * t31;
r(4) = m41 * t11 + m42 * t21 + m43 * t31;
r(5) = m11 * t12 + m12 * t22 + m13 * t32;
r(6) = m21 * t12 + m22 * t22 + m23 * t32;
r(7) = m31 * t12 + m32 * t22 + m33 * t32;
r(8) = m41 * t12 + m42 * t22 + m43 * t32;
r(9) = m11 * t13 + m12 * t23 + m13 * t33;
r(10) = m21 * t13 + m22 * t23 + m23 * t33;
r(11) = m31 * t13 + m32 * t23 + m33 * t33;
r(12) = m41 * t13 + m42 * t23 + m43 * t33;
if (r ~= m)
    r(13) = m(13);
    r(14) = m(14);
    r(15) = m(15);
    r(16) = m(16);
end
end
