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
        h = fspecial3.cubeAverage(siz);
              
    case 'sphereAverage'
        h = fspecial3.sphereAverage(siz);

    case 'gaussianFWHMhalfSize'
        h = fspecial3.gaussianFWHMhalfSize(siz);
             
    case 'laplacian'
        h = fspecial3.laplacian(siz);
     
    case 'log'
        h = fspecial3.log(siz);

    case 'dog'
        h = fspecial3.dog(siz);
                
     case 'gradient'    
        h = fspecial3.gradient(siz, dir);

    case 'gaussgradient1'
        h = fspecial3.gaussgradient1(siz, dir);
            
    case 'gaussgradient2'
        h = fspecial3.gaussgradient2(siz, dir);
     
    otherwise
        error('Unknown filter type.')
        
end

        