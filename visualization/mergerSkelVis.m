folder = '/home/mberning/Desktop/mergerNml/';

files = dir([folder '*.mat']);

close all;
figure('Position', [1921 1 1920 999]);
for i=1:length(files)
    subplot(5,2,2*i-1);
    load([folder files(i).name]);
    zeroIdx = ([results(:).segId1] == 0) & ([results(:).segId2]);
    temp = cellfun(@max, {results(~zeroIdx).prob}, 'UniformOutput', false);
    bar([sum(zeroIdx) sum(cellfun(@isempty, temp)) sum(~zeroIdx)-sum(cellfun(@isempty, temp))]);
    set(gca, 'XTick', [1 2 3]);
    set(gca, 'XTickLabel', {'zero hit' 'non-neighbour' 'normal edge'});
    ylabel('Fraction');
    title(strrep(strrep(files(i).name, ', ', ''), '_', ''));
    subplot(5,2,2*i);
    temp = cell2mat(temp(~cellfun(@isempty, temp)));
    hist(temp);
    xlabel('probability of normal edge');
end

export_fig(gcf, [folder 'visualization.svg']);
