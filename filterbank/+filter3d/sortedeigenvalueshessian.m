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
    Ixx = convn(I,   hxx{1}, 'same');
    Ixx = convn(Ixx, hxx{2}, 'same');
    Ixx = convn(Ixx, hxx{3}, 'same');

    Iyy = convn(I,   hyy{1}, 'same');
    Iyy = convn(Iyy, hyy{2}, 'same');
    Iyy = convn(Iyy, hyy{3}, 'same');

    Izz = convn(I,   hzz{1}, 'same');
    Izz = convn(Izz, hzz{2}, 'same');
    Izz = convn(Izz, hzz{3}, 'same');

    Ixy = convn(I,   hxy{1}, 'same');
    Ixy = convn(Ixy, hxy{2}, 'same');
    Ixy = convn(Ixy, hxy{3}, 'same');

    Ixz = convn(I,   hxz{1}, 'same');
    Ixz = convn(Ixz, hxz{2}, 'same');
    Ixz = convn(Ixz, hxz{3}, 'same');

    Iyz = convn(I,   hyz{1}, 'same');
    Iyz = convn(Iyz, hyz{2}, 'same');
    Iyz = convn(Iyz, hyz{3}, 'same');

    %Calculate eigenvectors and eigenvalues of the Hessian for each point
    Isize = size(I);
    imfeat = cell(3, 1);
    inputA = [Ixx(:)';Ixy(:)';Ixz(:)';
              Iyy(:)';Iyz(:)';Izz(:)';];
    roots = eig3S(inputA);
    Ieigen = roots';
    [~,sortIdx]=sort(abs(Ieigen),2);
    sortIdx = repmat((1:size(sortIdx,1))',1,3) + (sortIdx - 1)*size(sortIdx,1);
    Ieigen = Ieigen(sortIdx);
   
    imfeat{1} = reshape(Ieigen(:, 1), Isize);
    imfeat{2} = reshape(Ieigen(:, 2), Isize);
    imfeat{3} = reshape(Ieigen(:, 3), Isize);
end
