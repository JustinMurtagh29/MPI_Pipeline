function morphRcortex( param, nameAffMap, r )

    a = load([param.dataFolder param.affSubfolder nameAffMap '.mat'], 'class');
    r = param.r(r);
    if r ~= 0
        [x,y,z] = meshgrid(-r:r,-r:r,-r:r);
        se = (x/r).^2 + (y/r).^2 + (z/r).^2 <= 1;
        % Opening by reconstruction
        affEroded = imerode(a.class, se);
        affRecon = imreconstruct(affEroded, a.class);
        % Closing by reconstruction
        affReconDilated = imdilate(affRecon, se);
        affReconRecon = imreconstruct(imcomplement(affReconDilated), imcomplement(affRecon));
    else
        affReconRecon = imcomplement(a.class);
    end
    if ~exist([param.dataFolder param.outputSubfolder nameAffMap '/'], 'dir')
        mkdir([param.dataFolder param.outputSubfolder nameAffMap '/']);
    end
    parsave([param.dataFolder param.outputSubfolder filesep nameAffMap filesep 'MorphRecon' num2str(r) '.mat'], affReconRecon);



end

