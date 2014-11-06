function imfeat = gaussiansmoothedgradmagnitude(I, siz)

    hx=fspecial3.gaussgradient1(siz,1);
    hy=fspecial3.gaussgradient1(siz,2);
    hz=fspecial3.gaussgradient1(siz,3);

    Ix_x=convn(I,hx{1},'same');
    Ix_y=convn(Ix_x,hx{2},'same');
    Ix_z=convn(Ix_y,hx{3},'same');
    clear Ix_x Ix_y Ix_0
    Iy_x=convn(I,hy{1},'same');
    Iy_y=convn(Iy_x,hy{2},'same');
    Iy_z=convn(Iy_y,hy{3},'same');
    clear Iy_x Iy_y Iy_0
    Iz_x=convn(I,hz{1},'same');
    Iz_y=convn(Iz_x,hz{2},'same');
    Iz_z=convn(Iz_y,hz{3},'same');

    imfeat=sqrt(Ix_z.^2+Iy_z.^2+Iz_z.^2);
    clear Iz_x Iz_y Iz_0 Iz_z

end