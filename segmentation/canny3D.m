function final = canny3D(im, filsize, sigma, th_up, th_low)
% filsize = size filtermatrix (3x3x3 normally -> 3)
% sigma = 1,2,3...
% th_up ~ 0.07, th_low ~ 0.01
im = double(im);
%% 

hfil = floor(filsize/2);
[w,h,d] = size(im);
[x y z] = meshgrid(-hfil:hfil,-hfil:hfil,-hfil:hfil);

%Create Filter
fil_x = exp(-x.^2/(2*sigma^2))./(sigma*sqrt(2*pi)); clear x;
fil_y = exp(-y.^2/(2*sigma^2))./(sigma*sqrt(2*pi)); clear y;
fil_z = exp(-z.^2/(2*sigma^2))./(sigma*sqrt(2*pi)); clear z;
f = fil_x .* fil_y .* fil_z; clear fil_x; clear fil_y; clear fil_z;
f = f/sum(abs(f(:)));
imfil = imfilter(im,f,'replicate'); clear f;
[imfil_x, imfil_y, imfil_z] = gradient(imfil); clear imfil;

%%
%   Thinning (non-maximum suppression)
imfil_mag = sqrt(imfil_x.^2+imfil_y.^2+imfil_z.^2);

[w,h,d] = size(imfil_x);
[x,y,z] = meshgrid(1:h,1:w,1:d);

xi = x - imfil_x./imfil_mag;
yi = y - imfil_y./imfil_mag;
zi = z - imfil_z./imfil_mag;

imtemp = interp3(x,y,z,imfil_mag,xi,yi,zi);

xi = x + imfil_x./imfil_mag;
yi = y + imfil_y./imfil_mag;
zi = z + imfil_z./imfil_mag;

imtemp2 = interp3(x,y,z,imfil_mag,xi,yi,zi);

clear xi yi zi x y z imfil_x imfil_y imfil_z

%%
% hyseresis

im_sub = (imfil_mag > th_up);

[w,h,d] = size(im_sub);
index = find(im_sub);
old_index = [];

while length(old_index)~=length(index)
    old_index = index;
    x = mod(mod(index,w*h),h);
    y = ceil(mod(index,w*h)/h);
    z = ceil(index/w/h);
    for i = -1:1
        for j = -1:1
            for k = -1:1
                xtemp = x+i;
                ytemp = y+j;
                ztemp = z+k;
                ind = find(xtemp>=1 & xtemp<=w & ytemp>=1 & ytemp<=h & ztemp>=1 & ztemp<=d);
                xtemp = xtemp(ind);
                ytemp = ytemp(ind);
                ztemp = ztemp(ind);
                im(xtemp+(ytemp-1)*w+(ztemp-1)*h*w) = imfil_mag(xtemp+(ytemp-1)*w+(ztemp-1)*h*w)>th_low;
            end
        end
    end
    index = find(im);
end

im_sub = im_sub & (imtemp < imfil_mag)&(imcomplement(imtemp2 < imfil_mag));

clear im_th imfil_theta_z imtemp imtemp2 ind index xtemp ytemp ztemp 

final = im_sub;