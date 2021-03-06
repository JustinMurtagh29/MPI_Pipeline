function [x_shift1, y_shift1, x_shift2, y_shift2] = detectImageAlignmentIssues(p, planeIdx)

    raw1 = loadPlane(p.raw, p.bbox, planeIdx);
    raw2 = loadPlane(p.raw, p.bbox, planeIdx+1);
    xSteps = p.bbox(1,1)+80:40:p.bbox(1,2)-p.bbox(1,1)-80;
    ySteps = p.bbox(2,1)+80:40:p.bbox(2,2)-p.bbox(2,1)-80;
    [X,Y] = ndgrid(xSteps, ySteps);
    for i=1:length(X)
        for j=1:length(Y)
            fov1 = raw1(X(i)-80:X(i)+80,Y(j)-80:Y(j)+80);
            fov2 = raw2(X(i)-80:X(i)+80,Y(j)-80:Y(j)+80);
            % This first function still yields useless results
            %[x_shift1, y_shift1]= connectEM.xcorr2_fft(fov1, fov2);
            % This seems to work
            [x_shift1(i,j), y_shift1(i,j), x_shift2(i,j), y_shift2(i,j)]= connectEM.xcorr2_template(fov1, fov2);
       end
    end
end

function raw = loadPlane(raw, bbox, planeIdx)
    bbox(3,1:2) = [planeIdx planeIdx];
    raw = loadRawData(raw, bbox);
end

