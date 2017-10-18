function [x_shiftA, y_shiftA, x_shiftB, y_shiftB] = xcorr2_template(A, B)
    % Cross-correlation for 2 dimensional arrays
    % A and B should habe uneven number of elements in each dim and be square

    templateHalfSize = 50;
    A = single(A);
    B = single(B);

    % Do template matching in both directions to check for symetry
    midA = (size(A,1)-1)/2+1;
    midB = (size(B,1)-1)/2+1;
    % templates centered in the large image
    templateA = A(midA-templateHalfSize:midA+templateHalfSize,midA-templateHalfSize:midA+templateHalfSize);
    templateB = B(midB-templateHalfSize:midB+templateHalfSize,midB-templateHalfSize:midB+templateHalfSize);

    if numel(unique(templateA)) > 1
        xcorrB = normxcorr2(templateA, B);
        [x_max, y_max] = find(xcorrB == max(xcorrB(:)),1);
        % Minus to have shifts in the same direction
        x_shiftB = -(x_max-midB-templateHalfSize);
        y_shiftB = -(y_max-midB-templateHalfSize);
    else
        x_shiftB = 0;
        y_shiftB = 0;
    end

    if numel(unique(templateB)) > 1
        xcorrA = normxcorr2(templateB, A);
        [x_max, y_max] = find(xcorrA == max(xcorrA(:)),1);
        x_shiftA = x_max-midA-templateHalfSize;
        y_shiftA = y_max-midA-templateHalfSize;
    else
        x_shiftA = 0;
        y_shiftA = 0;
    end

end

