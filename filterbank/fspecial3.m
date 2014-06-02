function h = fspecial3(type,siz,dir)
% Created by Manuel Berning to produce 3D filters

%added: 
%dog: difference of gaussians, second gaussian has a width that is 2/3 of
%the first one
%gaussiangradient1/2: filter that applied to an image gives the
%first/second derivative wrt to dir (two direction to be specified for the
%second derivative) of the gaussian-smoothed image

if nargin < 3
    dir = [];
end


switch type
    case 'cubeAverage'
        
        h{1} = ones(1,siz,1)/prod([siz]);
        h{2} = ones(siz,1,1)/prod([siz]);
        h{3} = ones(1,1,siz)/prod([siz]);
              
    case 'sphereAverage'
        R = siz/2;
        R(R==0) = 1;
        h = ones(siz,siz,siz);
        siz = (siz-1)/2;
        [x,y,z] = ndgrid(-siz:siz,-siz:siz,-siz:siz);
        I = (x.*x/R^2+y.*y/R^2+z.*z/R^2)>1;
        h(I) = 0;
        h = h/sum(h(:));

    case 'gaussianFWHMhalfSize'
    
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
             h = {hx, hy, hz};
             
    case 'laplacian'
        if siz < 3
            error('Discrete Laplacian with size per dimension < 3 not defined');
        else
            siz = siz + [2 2 2];
            h = zeros(siz);
            h([1 end],:,:) = 1;
            h(:,[1 end],:) = 1;
            h(:,:,[1 end]) = 1;
            h = bwdist(h);
            h([1 end],:,:) = [];
            h(:,[1 end],:) = [];
            h(:,:,[1 end]) = [];
            if rem(siz,2) == 0  
                mid = (siz-2)/2;
                h(mid:mid+1, mid:mid+1, mid:mid+1) = 0;
                h(mid:mid+1, mid:mid+1, mid:mid+1) = repmat(-sum(h(:)),[2 2 2]);
                h = double(h);
            else
                mid = ceil((siz-2)/2);
                h(mid, mid, mid) = 0;
                h(mid, mid, mid) = sum(-h(:));
                h = double(h);
            end
        end
     
    case 'log'
        
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

    case 'dog'
        
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
                
     case 'gradient'    
        
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
        
    case 'gaussgradient1'
        
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
            
    case 'gaussgradient2'
        
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
     
    otherwise
        error('Unknown filter type.')
        
end

        