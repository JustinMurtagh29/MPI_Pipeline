function h = gradient(siz, dir)
if siz < 3
    error('Gradient Filter with size per dimension < 3 not defined');
elseif isempty(dir)
    error('Specifiy direction for gradient!');
else
    if rem(siz,2) == 0
        mid = siz/2-0.5;
        a= (-mid:mid)';
    else
        mid = floor(siz/2);
        a = (-mid:mid)';
    end
    b = ones(1,siz,1);
    switch dir
        case 1
            h{1} = b;
            h{2} = a;
            h{3} = ones(1,1,siz);
        case 2
            h{1} = b.';
            h{2} = a.';
            h{3} = ones(1,1,siz);
            % h = repmat(reshape(a * b, 1, siz, siz), [siz 1 1]);
        case 3
            h{1} = b.';
            h{2} = ones(1,1,siz);
            h{3} = a.';
        otherwise
            error('Only implemented for 3 axis aligned dimensions.')
    end
end
end

