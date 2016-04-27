function imfeat = sortedeigenvaluesstructure(I, siz, siz2)
    %calculate gradients with fspecial3 using size siz2
    hx=fspecial3.gradient(siz2,1);
    hy=fspecial3.gradient(siz2,2);
    hz=fspecial3.gradient(siz2,3);

    Ix_x=convn(I,hx{1},'same');
    Ix_y=convn(Ix_x,hx{2},'same');
    clear Ix_x
    Ix=convn(Ix_y,hx{3},'same');
    clear Ix_y
    Iy_x=convn(I,hy{1},'same');
    Iy_y=convn(Iy_x,hy{2},'same');
    clear Iy_x
    Iy=convn(Iy_y,hy{3},'same');
    clear Iy_y
    Iz_x=convn(I,hz{1},'same');
    Iz_y=convn(Iz_x,hz{2},'same');
    clear Iz_x
    Iz=convn(Iz_y,hz{3},'same');

    clear Iz_y hx hy hz

    %calculate elements of the structure tensor
    h=fspecial3.gaussianFWHMhalfSize(siz2); %window function

    Ixx1=convn(Ix.^2,h{1},'same');
    Ixx2=convn(Ixx1,h{2},'same');
    Ixx=convn(Ixx2,h{3},'same');
    clear Ixx1 Ixx2

    Ixy1=convn(Ix.*Iy,h{1},'same');
    Ixy2=convn(Ixy1,h{2},'same');
    Ixy=convn(Ixy2,h{3},'same');
    clear Ixy1 Ixy2

    Ixz1=convn(Ix.*Iz,h{1},'same');
    Ixz2=convn(Ixz1,h{2},'same');
    Ixz=convn(Ixz2,h{3},'same');
    clear Ixz1 Ixz2

    Iyy1=convn(Iy.^2,h{1},'same');
    Iyy2=convn(Iyy1,h{2},'same');
    Iyy=convn(Iyy2,h{3},'same');
    clear Iyy1 Iyy2

    Iyz1=convn(Iy.*Iz,h{1},'same');
    Iyz2=convn(Iyz1,h{2},'same');
    Iyz=convn(Iyz2,h{3},'same');
    clear Iyz1 Iyz2

    Izz1=convn(Iz.^2,h{1},'same');
    Izz2=convn(Izz1,h{2},'same');
    Izz=convn(Izz2,h{3},'same');               
    clear Izz1 Izz2 Iz Ix Iy

    %calculate the eigenvalues of the structure tensor
    [a,b,c]=size(I);
    imfeat=cell(3,1);
    newSize = [1 1 a*b*c];
    Ieigen = eig3([reshape(Ixx, newSize) reshape(Ixy, newSize) reshape(Ixz, newSize); ...
        reshape(Ixy, newSize) reshape(Iyy, newSize) reshape(Iyz, newSize); ...
        reshape(Ixz, newSize) reshape(Iyz, newSize) reshape(Izz, newSize)]);
    Ieigen = sort(Ieigen,2);
    imfeat{1}=reshape(Ieigen(:,1), [a b c]);
    imfeat{2}=reshape(Ieigen(:,2), [a b c]);
    imfeat{3}=reshape(Ieigen(:,3), [a b c]);  

end