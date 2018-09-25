function binIm = fillHoles(binIm, mask,oldmethod)
% this function takes mask and either 
% drill holes into 'outer hull' to be able to use imfill (oldmethod) or
% chooses a point most far away from any labeled object, fills from there
% and checks which parts cannot be reached ( = holes) which are then filled

binIm = or(binIm,mask);

if exist('oldmethod','var') && oldmethod
    %     holeLoc = round(size(binIm,1)/2);
    [~,holeLoc] = min(sum(binIm,2));
    
    if all(binIm(holeLoc,:))
        warning('Line where hole should be drilled is all black. No hole is drilled!')
    else
        idx = 1;
        while binIm(holeLoc,idx) == 1
            binIm(holeLoc,idx) = 0;
            idx = idx + 1;
        end
    end
    binIm = imfill(binIm, 'holes') & ~mask;
else
    
    dst = bwdist(binIm);
    [~,ind] = max(dst(:));
    [x,y] = ind2sub(size(dst),ind);
    binIm = imfill(binIm,[x,y])==0 | binIm & ~mask;
end
end