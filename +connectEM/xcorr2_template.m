function [x_shiftA, y_shiftA, x_shiftB, y_shiftB] = xcorr2_template(A, B)
    % Cross-correlation for 2 dimensional arrays
    % A and B should habe uneven number of elements in each dim and be square

    templateHalfSize = 50;
    A = single(A);
    B = single(B);

    % Do template matching in both directions to check for symetry
    midA = (size(A,1)-1)/2+1;
    midB = (size(B,1)-1)/2+1;
    % 11x11 templates centered in the large image
    templateA = A(midA-templateHalfSize:midA+templateHalfSize,midA-templateHalfSize:midA+templateHalfSize);
    templateB = B(midB-templateHalfSize:midB+templateHalfSize,midB-templateHalfSize:midB+templateHalfSize);
    xcorrA = normxcorr2(templateB, A);
    xcorrB = normxcorr2(templateA, B);

    [x_max, y_max] = find(xcorrA == max(xcorrA(:)));
    x_shiftA = x_max-midA-templateHalfSize;
    y_shiftA = y_max-midA-templateHalfSize;

    [x_max, y_max] = find(xcorrB == max(xcorrB(:)));
    % Minus to have shifts in the same direction
    x_shiftB = -(x_max-midB-templateHalfSize);
    y_shiftB = -(y_max-midB-templateHalfSize);

end

