function [x_max, y_max] = detectImageAlignmentIssues(p, planeIdx)

    raw1 = loadPlane(p.raw, p.bbox, planeIdx);
    raw2 = loadPlane(p.raw, p.bbox, planeIdx+1);
    xSteps = p.bbox(1,1)+80:40:p.bbox(1,2)-p.bbox(1,1)-80;
    ySteps = p.bbox(2,1)+80:40:p.bbox(2,2)-p.bbox(2,1)-80;
    [X,Y] = ndgrid(xSteps, ySteps);
    for i=1:length(X)
        for j=1:length(Y)
            fov1 = raw1(X(i)-80:X(i)+80,Y(j)-80:Y(j)+80);
            fov2 = raw2(X(i)-80:X(i)+80,Y(j)-80:Y(j)+80);
            xcorr = xcorr2_fft(fov1, fov2);
            [~, idx] = max(xcorr(:));
            [x_max(i,j), y_max(i,j)] = ind2sub(size(xcorr), idx);
        end
    end
end

function raw = loadPlane(raw, bbox, planeIdx)
    bbox(3,1:2) = [planeIdx planeIdx];
    raw = loadRawData(raw, bbox);
end

function result = xcorr2_fft(a, b)
    result = real(ifft2(fft2(a).*fft2(b(end:-1:1,end:-1:1))));
end

