function raw = copyRoi()
    % Copy part of 2012_09_28_ex145_07x2_new2 dataset availible on webKnossos to a new dataset to have the ROI for segmentation aligned to KNOSSOS cubes 

    % Define source of copy process
    chosenOne.root = '/gaba/u/kboerg/st07x2_new2/color/1/';
    chosenOne.prefix = '2012-09-28_ex145_07x2_new2_mag1'; 
    chosenOne.bbox = [1195, 1515, 115, 6690, 9945, 3420];
    % Convert to Matlab format (and add +1 offset between coordinate systems)
    chosenOne.bbox = Util.convertWebknossosToMatlabBbox(chosenOne.bbox);
    
    % Define destination of copy process
    destination.root = '/gaba/u/mberning/wkCubes/2012-09-28_ex145_07x2_ROI2016/color/1/';
    destination.prefix = '2012-09-28_ex145_07x2_ROI2016_mag1';
    destination.firstCoord = 129 - [25 25 10];
    
    % Uncomment for reading / writing .WKW files
    % chosenOne.backend = 'wkwrap';
    % destination.backend = 'wkwrap';
    
    if isfield(destination, 'backend') ...
            && strcmp(destination.backend, 'wkwrap')
        % initialize WKW dataset, if needed
        wkwInit('new', destination.root, 32, 32, 'uint8', 1);
    end

    % Read from one, write to wKcubes
    raw = loadRawData(chosenOne, chosenOne.bbox);
    saveRawData(destination, destination.firstCoord, raw);
end

