function imfeat = differenceofgaussians(I, siz)

    h = fspecial3.dog(siz);

    imfeatx = convn(I,h{1},'same');
    imfeaty = convn(imfeatx,h{2},'same');
    imfeat1 = convn(imfeaty,h{3},'same');

    imfeatx2 = convn(I,h{4},'same');
    imfeaty2 = convn(imfeatx2,h{5},'same');
    imfeat2 = convn(imfeaty2,h{6},'same');

    imfeat = imfeat1 - imfeat2;

end