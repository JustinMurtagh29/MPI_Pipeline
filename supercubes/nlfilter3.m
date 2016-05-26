function val = nlfilter3(data, fun, fS)

dS = size(data);
ri1=reshape(data,[dS(1) fS(2) dS(2)/fS(2) dS(3)]);
ri1=permute(ri1,[2 1 3 4]);
ri1=reshape(ri1,[fS(1)*fS(2) dS(1)*dS(2)/(fS(1)*fS(2)) fS(3) dS(3)/fS(3)]);
ri1=permute(ri1,[1 3 2 4]);
ri1=reshape(ri1,[fS(1)*fS(2)*fS(3) dS(1)*dS(2)*dS(3)/(fS(1)*fS(2)*fS(3))]);
ri1=fun(ri1);
val=reshape(ri1,[dS(1)/fS(1) dS(2)/fS(2) dS(3)/fS(3)]);

end

