function isosurfacesForPrinting(dendrites)

    load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat');
    p.local(end,:,:) = [];
    p.local(:,end,:) = [];
    p.local(:,:,end) = [];
    for i=1:length(dendrites)
        mask = Visualization.buildMaskLowRes(p, dendrites{i});
        Util.save(['/gaba/scratch/mberning/mask_' num2str(i, '%.2i') '.mat'], mask);
    end

end

