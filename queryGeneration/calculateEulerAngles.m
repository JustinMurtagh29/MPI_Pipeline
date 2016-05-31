function [phi, thetha, psi] = calculateEulerAngles(di)
    % Calculate angles (as deinfed in wK) in degrees from direction vector 

    % Make sure direction is normalized
    di = di ./ norm(di);
    % Get the LCS axis in GCS 
    [or1, or2] = findOrthogonals(di);
    % Calculate euler angles 
    [phi, thetha, psi] = LCS2Euler(di,or1,or2); 
    % Make angles positive
    phi = mod(phi+360,360);
    thetha = mod(thetha+360,360);
    psi = mod(psi+360,360);

end

function [orth1, orth2] = findOrthogonals(v)
    v = v ./ norm(v);
    if all(abs(v) == [0 0 1])
        orth1 = [1 0 -1].*v([3 2 1]);  
    else
        orth1 = [1 -1 0].*v([2 1 3]);
    end
    orth1 = orth1 ./ norm(orth1);
    orth2 = cross(v, orth1); 
end

function [phi, thetha, psi] = LCS2Euler(X1,Y1,Z1)
    Z1xy = sqrt(Z1(1).^2+Z1(2).^2);
    if (Z1xy > 1e-3)
        phi = atan2d(Y1(1)*Z1(2)-Y1(2)*Z1(1), X1(1)*Z1(2)-X1(2)*Z1(1));
        thetha = atan2d(Z1xy, Z1(3));
        psi = -atan2d(-Z1(1), Z1(2));
    else 
        phi = 0;
        if Z1(3) > 0
            thetha = 0;
        else
            thetha = 360;
        end
        psi = -atan2d(X1(2), X1(1));
    end
end
