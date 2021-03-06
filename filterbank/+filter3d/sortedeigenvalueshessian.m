function imfeat = sortedeigenvalueshessian(I, siz)
    %start by creating 6 filter for each entry in the Hessian

    %apply the filter to the image, ignoring boundaries, i.e. intensity outside
    %the image is set to zero
    
    hxx = fspecial3.gaussgradient2(siz, 11);
    Ixx = convn(I,   hxx{1}, 'same');
    Ixx = convn(Ixx, hxx{2}, 'same');
    Ixx = convn(Ixx, hxx{3}, 'same');
    clear hxx;

    hyy = fspecial3.gaussgradient2(siz, 22);
    Iyy = convn(I,   hyy{1}, 'same');
    Iyy = convn(Iyy, hyy{2}, 'same');
    Iyy = convn(Iyy, hyy{3}, 'same');
    clear hyy;

    hzz = fspecial3.gaussgradient2(siz, 33);
    Izz = convn(I,   hzz{1}, 'same');
    Izz = convn(Izz, hzz{2}, 'same');
    Izz = convn(Izz, hzz{3}, 'same');
    clear hzz;

    hxy = fspecial3.gaussgradient2(siz, 12);
    Ixy = convn(I,   hxy{1}, 'same');
    Ixy = convn(Ixy, hxy{2}, 'same');
    Ixy = convn(Ixy, hxy{3}, 'same');
    clear hxy;

    hxz = fspecial3.gaussgradient2(siz, 13);
    Ixz = convn(I,   hxz{1}, 'same');
    Ixz = convn(Ixz, hxz{2}, 'same');
    Ixz = convn(Ixz, hxz{3}, 'same');
    clear hxz;

    hyz = fspecial3.gaussgradient2(siz, 23);
    Iyz = convn(I,   hyz{1}, 'same');
    Iyz = convn(Iyz, hyz{2}, 'same');
    Iyz = convn(Iyz, hyz{3}, 'same');
    clear hyz;

    %Calculate eigenvectors and eigenvalues of the Hessian for each point
    Isize = size(I);
    imfeat = cell(3, 1);
    
    % prepare intput
    Ieigen = [Ixx(:)'; Ixy(:)'; Ixz(:)'; Iyy(:)'; Iyz(:)'; Izz(:)'];
    clear Ixx Ixy Ixz Iyy Iyz Izz;
    
    % compute eigenvalues
    Ieigen = eig3S(Ieigen)';
    Ieigen = single(Ieigen);
    
    % Eigen values are sorted based on their absolute values to be in consistency with classifiers trained on previously used 'eig' function
    [~, sortIds] = sort(abs(Ieigen), 2);
    
    % Each row of sortIdx matrix points to sorting order for each row of Ieigen. These indices are converted to linear indices to implement this on Ieigen
    linearIds = size(sortIds, 1);
    linearIds = bsxfun(@plus, (sortIds -1) * linearIds, (1:linearIds)');
    Ieigen = Ieigen(linearIds);
   
    imfeat{1} = reshape(Ieigen(:, 1), Isize);
    imfeat{2} = reshape(Ieigen(:, 2), Isize);
    imfeat{3} = reshape(Ieigen(:, 3), Isize);
end
