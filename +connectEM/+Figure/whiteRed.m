function cmap = whiteRed(n)
    % cmap = whiteRed(n)
    %   Generates a color map with `n` entries from white to red.
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    gamma = 2.2;
    
    cmap = zeros(n, 3);
    cmap(:, 1) = 1;
    
    alpha = linspace(0, 1, n);
    alpha = reshape(alpha, [], 1);
    
    cmap = ( ...
        (cmap .^ gamma) .* alpha ...
      + ([1, 1, 1] .^ gamma) .* (1 - alpha)...
      ) .^ (1 / gamma);
end
