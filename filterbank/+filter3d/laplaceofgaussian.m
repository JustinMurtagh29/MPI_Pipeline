function imfeat = laplaceofgaussian(I, siz)

    clear h  
    h=fspecial3.log(siz); 
    imfeat_x=convn(I,h{1}, 'same');
    imfeat_y=convn(imfeat_x,h{2},'same');
    imfeat=convn(imfeat_y,h{3},'same');

end