function sphere = makeSphere(size, r)
    % Based on code by Benedikt Staffler
    % From Util.getPointsInBall
    
    r2 = ceil(size);
    [xx, yy, zz] = meshgrid(-r2(1):r2(1), -r2(2):r2(2), -r2(3):r2(3));
    sphere = sqrt(xx .^ 2 + yy .^ 2 + (zz.*28./11.24) .^ 2) <= r;
end
