function imfeat = intensitygaussiansmoothed(I, siz)

    h=fspecial3.gaussianFWHMhalfSize(siz);             
    imfeatx = convn(I,h{1},'same');
    imfeaty = convn(imfeatx,h{2},'same');
    imfeat = convn(imfeaty,h{3},'same');
    clear imfeaty imfeatx

end