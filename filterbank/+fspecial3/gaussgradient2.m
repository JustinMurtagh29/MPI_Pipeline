function h = gaussgradient2(siz, dir)

if isempty(dir)
            error('Specify direction for gradient!');
        else 
            clear h;
            sig = siz/(4*sqrt(2*log(2)));
            size   = (siz-1)/2;
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
                    
            switch dir
                case 11
                   h{1}=-hx;
                   h{2}=(((y.').^2)-sig^2).*hy./sig^4;
                   h{3}=-hz;
                case 22
                   h{1}=((x.^2)-sig^2).*hx./sig^4;
                   h{2}=-hy;
                   h{3}=-hz;
                case 33
                   h{1}=-hx;
                   h{2}=-hy;
                   h{3}=((permute(z,[1,3,2]).^2)-sig^2).*hz./sig^4;
                case 12                                                     
                   h{1}=x.*hx./sig^2;
                   h{2}=-(y.').*hy./sig^2;
                   h{3}=-hz;
                case 13                                                    
                   h{1}=hx;
                   h{2}=-(x.').*hy./sig^2;
                   h{3}=-permute(z,[1,3,2]).*hz./sig^2;
                case 23                                                    
                   h{1}=-x.*hx./sig^2;
                   h{2}=hy;
                   h{3}=-permute(z,[1,3,2]).*hz./sig^2;
                otherwise
                    error('For a second derivative two directions have to be specified.')
            end
        end
end