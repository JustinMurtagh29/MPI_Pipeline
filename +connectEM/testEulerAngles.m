function result = testEulerAngles()

    tests = { [0 0 1; 0 0 NaN; 1 1 1], ...
        [0 0 -1; 180 0 NaN; 0 180 NaN], ...
        [0 1 0; 90, NaN, 0; 270 NaN 0], ...
        [0 -1 0; 270, NaN, 0; 270 NaN 180], ...
        [1 0 0; 90, NaN, 270; 270, NaN, 90], ...
        [-1 0 0; 90, NaN, 90; 270, NaN, 270]};

    for i=1:length(tests)
        [phi, thetha, psi] = calculateEulerAngles(tests{i}(1,:));
        thisTest = tests{i}(2:end,:);
        thisResult = [phi, thetha, psi];
        result(i,1:3)= thisResult;
        temp = thisTest'; 
        result(i,4:4+numel(thisTest)-1) = temp(:)';
        temp = bsxfun(@eq, thisTest, thisResult);
        temp(isnan(thisTest)) = true;
        result(i,10) = any(all(temp,2));
    end

end

function [phi, thetha, psi] = calculateEulerAngles(di)
    % Calculate angles (as deinfed in wK) in degrees from direction vector 

    % Make sure direction is normalized
    di = di ./ norm(di);
    % Get the LCS axis in GCS 
    %[or1, or2] = findOrthogonals(di);
    % Calculate euler angles
    %rotM = [or1;or2;di]'; 
    eulerAngles = diffToEulerAngle(di);
    %eulerAngles = rad2deg(rotm2eul(rotM'));
    phi = eulerAngles(1);
    thetha = eulerAngles(2);
    psi = eulerAngles(3);
    % Make angles positive for wK
    %phi = round(mod(phi+360,360));
    %thetha = round(mod(thetha+360,360));
    %psi = round(mod(psi+360,360));

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

