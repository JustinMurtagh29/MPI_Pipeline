function q = generateQueriesFromAgglos(p, segmentMeta, agglos, outputFolder)

    voxelSize = p.raw.voxelSize;
    % Decide which positions to querry and calculate some statistics
    [q.pos, q.dir, q.varianceExplained, q.voxelSize, q.lengthAlongPC1] = ...
        cellfun(@(x)determineQueryLocation(segmentMeta, x, voxelSize), agglos, 'uni', 0);
    % Sort out all queries based on some heuristics
    borderSize = round(1000 ./ voxelSize);
    bbox(:,1) = p.bbox(:,1) + borderSize';
    bbox(:,2) = p.bbox(:,2) - borderSize';
    q.outsideBBox = cellfun(@(x)any(bsxfun(@le, x, bbox(:,1)'),2) | any(bsxfun(@ge, x, bbox(:,2)'),2), q.pos, 'uni', 0);
    q.tooSmall = cellfun(@(x)x < 5000, q.voxelSize);
    q.tooShort = cellfun(@(x)x < 1500, q.lengthAlongPC1);
    q.tooUnstraight = cellfun(@(x)x(1) < .7, q.varianceExplained);
    q.exclude = arrayfun(@(x,y,z,t,u) x{:} | y | z | t, q.outsideBBox, q.tooShort, q.tooSmall, q.tooUnstraight, 'uni', 0);
    % 'Write' problems to flight mode webKNOSSOS
    extend = round(1000 ./ voxelSize);
    fid = fopen([outputFolder datestr(clock,30) '_flightTasks.txt'], 'w');
    for i=1:length(q.pos)
        for j=1:size(q.pos{i},1)
            [phi, theta, psi] = calculateEulerAngles(-q.dir{i}(j,:)); 
            q.angles{i}(j,:) = [phi theta psi];
            minPos = q.pos{i}(j,:) - extend;
            sizeBbox = 2*extend;
            taskString = ['2012-09-28_ex145_07x2_ROI2017,56d6a7c6140000d81030701e,focus_flight,1,' ...
                num2str(q.pos{i}(j,1)-1) ',' num2str(q.pos{i}(j,2)-1) ',' num2str(q.pos{i}(j,3)-1) ',' ...
                num2str(phi) ',' num2str(theta) ',' num2str(psi) ',1,Tracing crew,' ...
                num2str(minPos(1)) ',' num2str(minPos(2)) ',' num2str(minPos(3)) ',' ...
                num2str(sizeBbox(1)) ',' num2str(sizeBbox(2)) ',' num2str(sizeBbox(3)) ',' 'L4_focus_flight_2017_test'];
            if ~q.exclude{i}(j)
                fprintf(fid, '%s\n', taskString);
            end
        end
    end
    fclose(fid);

end

function [p, d, vE, s, l] = determineQueryLocation(segmentMeta, agglo, voxelSize)
    % Determine query position, direction and edge based on PCA

    thesePoints = segmentMeta.point(:,agglo)';
    [coeff,score,latent] = pca(bsxfun(@times, thesePoints, voxelSize));
    [~,minIdx] = min(score(:,1));
    [~,maxIdx] = max(score(:,1)); 
    p(1,:) = thesePoints(minIdx,:);
    p(2,:) = thesePoints(maxIdx,:);
    d(1,:) = -coeff(:,1);
    d(2,:) = coeff(:,1);
    vE = latent./sum(latent);
    % Determine number of voxel for current agglo
    s = sum(segmentMeta.voxelCount(agglo));
    % Determine length of agglomeration along PC1
    l = sqrt(sum(((p(2,:) - p(1,:)).*voxelSize).^2));

end

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

