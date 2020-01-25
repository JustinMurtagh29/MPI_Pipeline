function cmap = redBlue(n)
    % cmap = redBlue(n)
    %   Generates a color map with `n` entries from blue over white to red.
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    gamma = 2.2;
    
    cmap = zeros(n, 3);
    cmap(1:ceil(n / 2), 3) = 1;
    cmap(ceil(n / 2 + 1):end, 1) = 1;
    
    alpha = [ ...
        linspace(1, 0, ceil(n / 2)), ...
        linspace(0, 1, n - ceil(n / 2))];
    alpha = reshape(alpha, [], 1);
    
    cmap = ( ...
        (cmap .^ gamma) .* alpha ...
      + ([1, 1, 1] .^ gamma) .* (1 - alpha)...
      ) .^ (1 / gamma);
end
