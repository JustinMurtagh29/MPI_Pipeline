function Pb = Pbfilter(im)

[f_e_x,f_e_y,f_e_z] = gradient(gradient((exp(-im.^2./2))));
f_even = (f_e_x + f_e_y + f_e_z)./3;
clear f_e_x f_e_y f_e_z

f_odd = zeros(size(im));
for i = 1 : size(im,3)
    f_odd(:,:,i) = imag(hilbert(im(:,:,i)));
end   

Pb = ((im.*f_even).^2 + (im.*f_odd).^2);

end