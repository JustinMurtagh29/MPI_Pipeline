function raw = copyRoi()
    % Copy part of 2012_09_28_ex145_07x2_new2 dataset availible on webKnossos to a new dataset to have the ROI for segmentation aligned to KNOSSOS cubes 

    % Define source of copy process
    chosenOne.root = '/gaba/u/alik/wKcubes/2017-04-03_ppcAK99_76x79/color/1/';
    chosenOne.prefix = '2017-04-03_ppcAK99_76x79_mag1'; 
    chosenOne.bbox = [1075, 1330, 0, 5880, 7585, 4704];
    % Convert to Matlab format (and add +1 offset between coordinate systems)
    chosenOne.bbox = Util.convertWebknossosToMatlabBbox(chosenOne.bbox);
    
    % Define destination of copy process
    destination.root = '/gaba/u/alik/wKcubes/2017-04-03_ppcAK99_76x79_preprocessing/color/1/';
    destination.prefix = '2017-04-03_ppcAK99_76x79_preprocessing_mag1';
    destination.firstCoord = 129 - [25 25 10];

    % Read from one, write to wKcubes
    raw = readKnossosRoi(chosenOne.root, chosenOne.prefix, chosenOne.bbox);
    writeKnossosRoi(destination.root, destination.prefix, destination.firstCoord, raw);

end

