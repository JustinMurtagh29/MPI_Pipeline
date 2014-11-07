function h = log(siz)

    if siz < 3
        error('Laplacian of Gaussian Filter with size per dimension < 3 not defined');
    else
        sig = siz/(4*sqrt(2*log(2)));
        size   = (siz-1)/2;
        x = -size:size;
        y = -size:size;
        z = -size:size;
        hx = exp(-(x.*x/2/sig^2));
        hy = exp(-(y.*y/2/sig^2)).';
        hz = ones(1,1,5);
        for i = 1:siz
            hz(1,1,i) = exp(-(z(i)^2/2/sig^2));
        end             
        arg = (x.*x/sig^4 + y.*y/sig^4 + z.*z/sig^4 - ...
            (1/sig^2 + 1/sig^2 + 1/sig^2));
        hx = arg.*hx;
        hx = hx/sum(hx(:));
        hy = hy/sum(hy(:));
        hz = hz/sum(hz(:));
        h = {hx, hy,hz};            
        
    end
end