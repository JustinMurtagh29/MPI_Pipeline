%folder = '/gaba/u/mberning/results/pipeline/20141007T094904/local/';
%foldersToLink = dir([folder 'x*']);

%for i=1:length(foldersToLink)
%    x = str2double(foldersToLink(i).name(2:5));
%    if x > 1 && x < 15
%        system(['ln -s ' folder foldersToLink(i).name ' /gaba/u/mberning/results/pipeline/20151111T183414/local/' ...
%            strrep(foldersToLink(i).name, ['x' num2str(x, '%.4i')], ['x' num2str(x-1, '%.4i')])]);
%    end
%end

folder = '/gaba/u/mberning/results/pipeline/20141007T094904/correspondences/';
filesToLink = dir([folder '*global.mat']);

for i=1:length(filesToLink)
    x1 = str2double(filesToLink(i).name(1:2));
    x2 = str2double(filesToLink(i).name(7:8));
    if x1 > 1 && x1 < 15 && x2 > 1 && x2 < 15
        system(['ln -s ' folder filesToLink(i).name ' /gaba/u/mberning/results/pipeline/20151111T183414/correspondences/' ...
            num2str(x1-1, '%.2i') filesToLink(i).name(3:6) num2str(x2-1, '%.2i') filesToLink(i).name(9:end)]);

    end
end

