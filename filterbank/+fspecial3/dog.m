function h = dog(siz)
    sig1=siz/(4*sqrt(2*log(2)));
    sig2=sig1*0.66;
    dsig=sig1-sig2;

    size   = (siz-1)/2;
    x = -size:size;
    y = -size:size;
    z = -size:size;

    hx = exp(-(x.*x/2/sig1^2));
    hx = hx/sum(hx(:));
    hy = exp(-(y.*y/2/sig1^2)).';
    hy = hy/sum(hy(:));
    hz = exp(-(z.^2/2/sig1^2));
    hz = permute(hz, [1,3,2]);
    hz = hz/sum(hz(:));

    hx2 = exp(-(x.*x/2/sig2^2));
    hx2 = hx2/sum(hx2(:));
    hy2 = exp(-(y.*y/2/sig2^2)).';
    hy2 = hy2/sum(hy2(:));
    hz2 = exp(-(z.*z/2/sig2^2));
    hz2 = permute(hz2, [1,3,2]);
    hz2 = hz2/sum(hz2(:));

    h{1}=hx;
    h{2}=hy;
    h{3}=hz;
    h{4}=hx2;
    h{5}=hy2;
    h{6}=hz2;

end