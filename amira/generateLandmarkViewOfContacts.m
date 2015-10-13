addpath(genpath('/home/mberning/code/SegEM/auxiliaryMethods/'));

folder = '/home/mberning/Desktop/20150902T192556/';
files = dir([folder '*.nml']);

for i=1:length(files)  
    skel = parseNml([folder files(i).name]);
    nrContacts = 0;
    clear coords;
    for j=1:length(skel)
       if strncmp(skel{j}.name,'Contact ID', 10)
           nrContacts = nrContacts + 1;
           coords(nrContacts,:) = mean(skel{j}.nodes(:,1:3));
       end
    end
    % Write to file
    if exist('coords', 'var')
        fid = fopen([folder strrep(files(i).name, '.nml', '.landmark')], 'w+');
        fprintf(fid, '%s\n\n\n', '# AmiraMesh 3D ASCII 2.0');
        fprintf(fid, '%s\n\n', ['define Markers ' num2str(nrContacts)]);
        fprintf(fid, '%s\n\t%s\n\t%s\n\n', 'Parameters {', 'NumSets 1,', 'ContentType "LandmarkSet"}');
        fprintf(fid, '%s\n\n', 'Markers { float[3] Coordinates } @1');
        fprintf(fid, '%s\n', '# Data section follows\n');
        fprintf(fid, '%s\n', '@1');
        for j=1:size(coords,1)
            fprintf(fid, '%21.15e %21.15e %21.15e \n', coords(j,1).*11.24, coords(j,2).*11.24, coords(j,3).*28);
        end
        fclose(fid); 
    end
end
