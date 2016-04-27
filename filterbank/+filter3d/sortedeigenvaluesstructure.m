function imfeat = sortedeigenvaluesstructure(I, siz, siz2)
    % calculate gradients with fspecial3 using size siz2
    hx =fspecial3.gradient(siz2, 1);
    hy =fspecial3.gradient(siz2, 2);
    hz =fspecial3.gradient(siz2, 3);

    Ix = convn(I,  hx{1}, 'same');
    Ix = convn(Ix, hx{2}, 'same');
    Ix = convn(Ix, hx{3}, 'same');
    
    Iy = convn(I,  hy{1}, 'same');
    Iy = convn(Iy, hy{2}, 'same');
    Iy = convn(Iy, hy{3}, 'same');
    
    Iz = convn(I,  hz{1}, 'same');
    Iz = convn(Iz, hz{2}, 'same');
    Iz = convn(Iz, hz{3}, 'same');

    clear Iz_y hx hy hz

    % calculate elements of the structure tensor
    % window function
    h = fspecial3.gaussianFWHMhalfSize(siz2);

    Ixx = convn(Ix .^ 2,  h{1}, 'same');
    Ixx = convn(Ixx,      h{2}, 'same');
    Ixx = convn(Ixx,      h{3}, 'same');

    Ixy = convn(Ix .* Iy, h{1}, 'same');
    Ixy = convn(Ixy,      h{2}, 'same');
    Ixy = convn(Ixy,      h{3}, 'same');

    Ixz = convn(Ix .* Iz, h{1}, 'same');
    Ixz = convn(Ixz,      h{2}, 'same');
    Ixz = convn(Ixz,      h{3}, 'same');

    Iyy = convn(Iy .^ 2,  h{1}, 'same');
    Iyy = convn(Iyy,      h{2}, 'same');
    Iyy = convn(Iyy,      h{3}, 'same');

    Iyz = convn(Iy .* Iz, h{1}, 'same');
    Iyz = convn(Iyz,      h{2}, 'same');
    Iyz = convn(Iyz,      h{3}, 'same');

    Izz = convn(Iz .^ 2,  h{1}, 'same');
    Izz = convn(Izz,      h{2}, 'same');
    Izz = convn(Izz,      h{3}, 'same');

    % calculate the eigenvalues of the structure tensor
    Isize = size(I);
    imfeat = cell(3, 1);
    Ieigen = eig3S([Ixx(:)'; Ixy(:)'; Ixz(:)'; Iyy(:)'; Iyz(:)'; Izz(:)'])';
    clear Ixx Ixy Ixz Iyy Iyz Izz
    % Eigen values are sorted based on their absolute values to be in consistency with classifiers trained on previously used 'eig' function
    [~,sortIds]=sort(abs(Ieigen),2);
    % Each row of sortIdx matrix points to sorting order for each row of Ieigen. These indices are converted to linear indices to implement this on Ieigen
    sortRows = size(sortIds,1);
    linearIds = bsxfun(@plus, (sortIds -1)*sortRows,(1:sortRows)');
    Ieigen = Ieigen(linearIds);

    imfeat{1} = reshape(Ieigen(:, 1), Isize);
    imfeat{2} = reshape(Ieigen(:, 2), Isize);
    imfeat{3} = reshape(Ieigen(:, 3), Isize);
end
