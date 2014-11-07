function h = sphereAverage(siz)
    R = siz/2;
    R(R==0) = 1;
    h = ones(siz,siz,siz);
    siz = (siz-1)/2;
    [x,y,z] = ndgrid(-siz:siz,-siz:siz,-siz:siz);
    I = (x.*x/R^2+y.*y/R^2+z.*z/R^2)>1;
    h(I) = 0;
    h = h/sum(h(:));
end