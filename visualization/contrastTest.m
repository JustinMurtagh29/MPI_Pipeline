function contrastTest(raw, fieldOfView)

MD = imfilter(raw,fspecial('average',fieldOfView),'same');
SDD = stdfilt(raw, ones(fieldOfView));

figure('Position', [1 41 1600 784]); 
imagesc(raw);
axis off; axis equal;
colormap('gray');
title('Original Image');
figure('Position', [1 41 1600 784]); 
imagesc(MD);
axis off; axis equal;
colormap('gray');
title('Mean Image');
figure('Position', [1 41 1600 784]); 
imagesc(SDD);
axis off; axis equal;
colormap('gray');
title('STD Image');

end
