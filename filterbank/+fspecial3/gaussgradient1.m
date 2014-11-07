function h = gaussgradient1(siz, dir)

if isempty(dir)
            error('Specify direction for gradient!');
        else 
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
                case 1
                    h{1}=-hx;
                    h{2}=-(x.').*hy./sig^2;
                    h{3}=-hz;
                case 2
                    h{1}=-y.*hx./sig^2;            
                    h{2}=-hy;
                    h{3}=-hz;
                case 3                             
                    h{1}=-hx;
                    h{2}=-hy;
                    h{3}=-(permute(z,[1,3,2])).*hz./sig^2;
                otherwise
                    error('Only implemented for 3 axis aligned dimensions.')
                    
            end
        end

end