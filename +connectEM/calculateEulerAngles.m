function [phi, thetha, psi] = calculateEulerAngles(di, voxelSize)
    % Calculate angles (as deinfed in wK) in degrees from direction vector

    % Make sure direction is normalized
    di = di .* voxelSize;
    di = di ./ norm(di);
    % Calculate euler angles, from wk-paper -> RESCOPaaS
    eulerAngles = connectEM.diffToEulerAngle(di);
    phi = round(eulerAngles(1));
    thetha = round(eulerAngles(2));
    psi = round(eulerAngles(3));

end

