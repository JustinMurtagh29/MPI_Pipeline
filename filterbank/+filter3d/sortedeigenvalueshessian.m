function imfeat = sortedeigenvalueshessian(I, siz)
    %start by creating 6 filter for each entry in the Hessian
    hxx=fspecial3.gaussgradient2(siz,11);
    hyy=fspecial3.gaussgradient2(siz,22);
    hzz=fspecial3.gaussgradient2(siz,33);
    hxy=fspecial3.gaussgradient2(siz,12);
    hxz=fspecial3.gaussgradient2(siz,13);
    hyz=fspecial3.gaussgradient2(siz,23);

    %apply the filter to the image, ignoring boundaries, i.e. intensity outside
    %the image is set to zero
    Ixx_x=convn(I,hxx{1},'same');
    Ixx_y=convn(Ixx_x,hxx{2},'same');
    Ixx=convn(Ixx_y,hxx{3},'same');

    Iyy_x=convn(I,hyy{1},'same');
    Iyy_y=convn(Iyy_x,hyy{2},'same');
    Iyy=convn(Iyy_y,hyy{3},'same');
    clear Iyy_x Iyy_y Ixx_x Ixx_y

    Izz_x=convn(I,hzz{1},'same');
    Izz_y=convn(Izz_x,hzz{2},'same');
    Izz=convn(Izz_y,hzz{3},'same');

    Ixy_x=convn(I,hxy{1},'same');
    Ixy_y=convn(Ixy_x,hxy{2},'same');
    Ixy=convn(Ixy_y,hxy{3},'same');
    clear Izz_x Izz_y Ixy_x Ixy_y

    Ixz_x=convn(I,hxz{1},'same');
    Ixz_y=convn(Ixz_x,hxz{2},'same');
    Ixz=convn(Ixz_y,hxz{3},'same');

    Iyz_x=convn(I,hyz{1},'same');
    Iyz_y=convn(Iyz_x,hyz{2},'same');
    Iyz=convn(Iyz_y,hyz{3},'same');
    clear Ixz_x Ixz_y Iyz_x Iyz_y

    %Calculate eigenvectors and eigenvalues of the Hessian for each point
    [a,b,c]=size(I);
    imfeat=cell(3,1);
    newSize = [1 1 a*b*c];
    Ieigen= eig3([reshape(Ixx, newSize) reshape(Ixy, newSize) reshape(Ixz, newSize); ...
        reshape(Ixy, newSize) reshape(Iyy, newSize) reshape(Iyz, newSize); ...
        reshape(Ixz, newSize) reshape(Iyz, newSize) reshape(Izz, newSize)]);
    Ieigen = sort(Ieigen,2);
    imfeat{1}=reshape(Ieigen(:,1), [a b c]);
    imfeat{2}=reshape(Ieigen(:,2), [a b c]);
    imfeat{3}=reshape(Ieigen(:,3), [a b c]); 

end