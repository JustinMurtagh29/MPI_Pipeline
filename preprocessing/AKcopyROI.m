function [] = AKcopyRoi()
    % Copy part of 2012_09_28_ex145_07x2_new2 dataset availible on webKnossos to a new dataset to have the ROI for segmentation aligned to KNOSSOS cubes 

    % Define source of copy process
    chosenOne.root = '/gaba/u/alik/wKcubes/2017-04-03_ppcAK99_76x79/color/1/';
    chosenOne.prefix = '2017-04-03_ppcAK99_76x79_mag1'; 
    chosenOne.bbox = [1075, 1330, 0, 5880, 7585, 4704];
    % Convert to Matlab format (and add +1 offset between coordinate systems)
    chosenOne.bbox = Util.convertWebknossosToMatlabBbox(chosenOne.bbox);
    
    % Define destination of copy process
    destination.root = '/gaba/u/alik/wKcubes/2017-04-03_ppcAK99_76x79_ROI/color/1/';
    destination.prefix = '2017-04-03_ppcAK99_76x79_ROI_mag1';
    zCoords = chosenOne.bbox(3,1):chosenOne.bbox(3,2);
    lastOfKnossosCube = [find(mod(zCoords,512) == 0) length(zCoords)];
    numberValidInCube = [lastOfKnossosCube(1) diff(lastOfKnossosCube)];
    zCoords = mat2cell(zCoords, 1, numberValidInCube);
    zCoords(numberValidInCube == 0) = [];
    thisSliceBbox=chosenOne.bbox; 
     for  i=4:length(zCoords)
        disp(num2str(i))
        thisSliceBbox(3,:) = [zCoords{i}(1) zCoords{i}(end)];
        disp('starting to read')
        % Read from one, write to wKcubes
        raw = readKnossosRoi(chosenOne.root, chosenOne.prefix, thisSliceBbox);
        disp('finished reading, starting to write')
        destination.firstCoord = 129 - [25 25 10];
        destination.firstCoord(1,3) =destination.firstCoord(1,3)+zCoords{i}(1)-1;
        writeKnossosRoi(destination.root, destination.prefix, destination.firstCoord, raw);
    end
end

