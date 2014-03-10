%% set all markers to a certain size (here 10)
set(findall(gcf,'-property','MarkerSize'),'MarkerSize',10);
% set all text to a certain size (here 20)
set(findall(gcf,'-property','FontSize'),'FontSize',20);
set(findall(gcf,'-property','LineWidth'),'LineWidth',2);
grid off; box off;

%% Set y axis limits
% ylim([1e2 1e5])
% 
% % Set x axis limits
% xlim([1e3 1e7])

%% Label change from nm to microns
set(gca, 'XTickLabel', {'1' '10' '100' '1000' '10000'});
set(gca, 'YTickLabel', {'0.1' '1' '10' '100'});
xlabel('average distance between merger [micron]');
ylabel('average distance between splits [micron]');

%% Some titles to choose from/add to
title('Split and merger rates cortex before NR');
title('Split and merger rates retina before NR');
title('Split and merger rates cortex after NR');
title('Split and merger rates retina after NR');
title('Comparison of test performance after NR');
title('Sampling test series (color indicates node density)');

%%
h = get(gcf, 'Children');
ax = get(h, 'Children');    
for i=1:1116
    set(ax(i), 'MarkerEdgeColor', [0 .6 0]);
end
for i=1117:2232
    set(ax(i), 'MarkerEdgeColor', [0 0 0]);
end