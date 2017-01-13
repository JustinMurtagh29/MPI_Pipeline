function smoothed = smooth3Aniso( raw, size, sig )
% Smooth 3 from MATLAB adaptation to allow for Gaussian sigma dependent on
% dimension (Matlab just allows setting filter size but NOT sigma per
% dimension (which seems weird). This only implements 'gaussian'
% size & sig should be 3 element vectors and size contain odd integers

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
    smoothed = convn(raw,hx,'same');
    smoothed = convn(smoothed,hy,'same');
    smoothed = convn(smoothed,hz,'same');
    
end

