function cmap = redBlue(n)
    % cmap = redBlue(n)
    %   Generates a color map with `n` entries from blue over white to red.
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    m = round(0.5 * n);
    cmap = nan(n, 3);
    
    cmap(:, 1) = min((1:n) / m, 1);
    cmap(:, 2) = 1 - abs(1 - (1:n) / m);
    cmap(:, 3) = min(2 - (1:n) / m, 1);
end
