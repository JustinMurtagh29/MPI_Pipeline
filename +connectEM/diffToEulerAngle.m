function eulerA = diffToEulerAngle(diffM)
% Kevin M. Boergens, 2016

    rotationMatrix = connectEM.angle_fudge2(diffM);
    rotationMatrixT = reshape((reshape(rotationMatrix,4,4)'),1, 16);
    eulerA = mod(connectEM.my2E(rotationMatrixT) * 180 / pi + 360, 360);

end

