function h = gaussianFWHMhalfSize(siz)
    sig = siz/(4*sqrt(2*log(2)));
    size = (siz-1)/2;
    x = -size:size;
    y = -size:size;
    z = -size:size;

    hx = exp(-(x.*x/2/sig^2));
    hx = hx/sum(hx(:));
    hy = exp(-(y.*y/2/sig^2)).';
    hy = hy/sum(hy(:));
    hz = exp(-(z.^2/2/sig^2));
    hz = permute(hz, [1,3,2]);
    hz = hz/sum(hz(:));
    h = {hx, hy, hz};
end