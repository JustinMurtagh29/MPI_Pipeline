function gradient = gradient3Aniso( raw, size, sig )
% Calculates gradient magnitude for anisotropic data

    % Make grid with 0 in the middle
    size = (size-1)/2;
    x = -size(1):size(1);
    y = -size(2):size(2);
    z = -size(3):size(3);

    % Smoothing for each direction
    hx = exp(-(x.*x/2/sig(1)^2));
    hx = hx/sum(hx(:));
    hy = exp(-(y.*y/2/sig(2)^2)).';
    hy = hy/sum(hy(:));
    hz = exp(-(z.^2/2/sig(3)^2));
    hz = permute(hz, [1,3,2]);
    hz = hz/sum(hz(:));
    
    % Apply filter
    gradient_x = convn(raw,x.*hx,'same');
    gradient_y = convn(raw,y'.*hy,'same');
    gradient_z = convn(raw,permute(z, [1,3,2]).*hz,'same');
    gradient = sqrt(gradient_x.^2+gradient_y.^2+gradient_z.^2);
end

