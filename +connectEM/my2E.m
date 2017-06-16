function euler = my2E(rotation_matrix)
% generate euler angles from rotation matrix (intermediate as quaternion)
% adapted from THREEjs setFromQuaternion and setFromRotationMatrix
% Kevin M. Boergens, 2016
sx = norm(rotation_matrix(1 : 3));
sy = norm(rotation_matrix(5 : 7));
sz = norm(rotation_matrix(9 : 11));

if (det(reshape(rotation_matrix, 4, 4)) < 0)
    sx = -sx;
end


invSX = 1 / sx;
invSY = 1 / sy;
invSZ = 1 / sz;

rotation_matrix(1 : 3) = rotation_matrix(1 : 3) * invSX;
rotation_matrix(5 : 7) = rotation_matrix(5 : 7) * invSY;
rotation_matrix(9 : 11) = rotation_matrix(9 : 11) * invSZ;

m11 = rotation_matrix(1);
m12 = rotation_matrix(5);
m13 = rotation_matrix(9);
m21 = rotation_matrix(2);
m22 = rotation_matrix(6);
m23 = rotation_matrix(10);
m31 = rotation_matrix(3);
m32 = rotation_matrix(7);
m33 = rotation_matrix(11);


trace = m11 + m22 + m33;

if (trace > 0)
    
    s = 0.5 / sqrt(trace + 1.0);
    
    w = 0.25 / s;
    x = (m32 - m23) * s;
    y = (m13 - m31) * s;
    z = (m21 - m12) * s;
    
else if (m11 > m22 && m11 > m33)
        
        s = 2.0 * sqrt(1.0 + m11 - m22 - m33);
        
        w = (m32 - m23) / s;
        x = 0.25 * s;
        y = (m12 + m21) / s;
        z = (m13 + m31) / s;
        
    else if (m22 > m33)
            
            s = 2.0 * sqrt(1.0 + m22 - m11 - m33);
            
            w = (m13 - m31) / s;
            x = (m12 + m21) / s;
            y = 0.25 * s;
            z = (m23 + m32) / s;
            
        else
            s = 2.0 * sqrt(1.0 + m33 - m11 - m22);
            
            w = (m21 - m12) / s;
            x = (m13 + m31) / s;
            y = (m23 + m32) / s;
            z = 0.25 * s;
            
        end
    end
end
sqx = x * x;
sqy = y * y;
sqz = z * z;
sqw = w * w;
clamp = @(v)max(min(v, 1), -1);

euler = [atan2(2 * (x * w - y * z), (sqw - sqx - sqy + sqz))
    asin(clamp(2 * (x * z + y * w)))
    atan2(2 * (z * w - x * y), (sqw + sqx - sqy - sqz)) - pi];
