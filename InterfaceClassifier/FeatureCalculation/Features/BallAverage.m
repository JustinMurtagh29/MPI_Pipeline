function imfeat = BallAverage( raw, radius )
%BALLAVERAGE Apply a ball average filter to 3D volume (image stack).

[x,y,z] = meshgrid(-radius:radius,-radius:radius,-ceil(radius/2):ceil(radius/2));
h = x.^2 + y.^2 + (2.*z).^2 <= (radius)^2;
h = h/sum(h(:));
imfeat = imfilter(raw,double(h));


end

