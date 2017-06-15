folder = {'D:\gaba02\gaba\scratch\kboerg\eval_agglo\aggloPostHoc2\axons1\', 'D:\gaba02\gaba\scratch\kboerg\eval_agglo\rungrid5_564\axons1\'}
titles = {'posthoc', 'original'}
%load('D:\gaba02\gaba\scratch\kboerg\tempaadd.mat');
close all
for idx = 1 : 10
    figure(idx);
    for idx3 = 1 : 2
        liste = dir([folder{idx3} num2str(idx) '_*.nml']);
    
        subplot(1,2,idx3)
        title(titles{idx3});
        hold on;
        %tempff{idx} = setdiff(tempdd{idx}, tempaa{idx}, 'rows');
        %scatter3(tempaa{idx}(:,1)*11.24/1000, tempaa{idx}(:,2)*11.24/1000, tempaa{idx}(:,3)*28/1000, 'x');
        %scatter3(tempff{idx}(:,1)*11.24/1000, tempff{idx}(:,2)*11.24/1000, tempff{idx}(:,3)*28/1000, '+');
        %legend('axon', 'single');
        colors = {'r','g','b','y','c','m'};
        for idx2 = 1 : length(liste)
            skel = skeleton([folder{idx3} liste(idx2).name]);
            xx = skel.nodes{1}(:, 1);
            yy = skel.nodes{1}(:, 2);
            zz = skel.nodes{1}(:, 3);
            plot3(11.24/1000*xx(skel.edges{1})',11.24*yy(skel.edges{1})'/1000,28/1000*zz(skel.edges{1})', colors{mod(idx2 - 1, 6) + 1}, 'LineWidth', 2)
            %plotSkelTraceNew(skel.reverse(), colors{mod(idx2 - 1, 6) + 1}, 2);
        end
         skel2= skeleton([folder{idx3} 'skel_' num2str(idx) '.nml']);
    plotSkelTraceNew(skel2.reverse(),'k', 1);
   
    end
    
end
folder2 = 'D:\tempmay\10skels\';
mkdir(folder2)
for idx = 1 : 10
    f=figure(idx)
    f.Renderer='Painters';
    set(f,'PaperOrientation','landscape');
    %set(f,'PaperUnits','centimeter');
    %set(f, 'PaperSize',[200,200]);
    %set(f, 'PaperPositionMode','manual');
    %set(f, 'PaperPosition',[10,10,180,180]);
    
    subplot(1,2,1)
    daspect([1,1,1])
    view(3)
    subplot(1,2,2)
    daspect([1,1,1])
    view(3)
    saveas(gcf, [folder2 'overlayPostHoc' num2str(idx) '.pdf'])
    %saveas(gcf, [folder2 'overlay' num2str(idx) '.fig'])
end
%%
for idx = 1 : 10
    figure(idx)
    scatter(1 -VV{idx}, VVrank{idx}+1)
    xlabel('1 - prob')
    ylabel('rank')
    set(gca, 'XScale','log', 'YScale', 'log');
    title(['axon ' num2str(idx)])
    saveas(gcf, [folder2 'rank' num2str(idx) '.pdf'])
end
%%
todo=cell(1,10);
for idx = 1 : 10
    for idx2 = 1: length(ytrue.VV{idx})
        reduction= length(intersect(ytrue.VVrankInfo{idx}{idx2}, ytrue.VVinfo{idx}));
        todo{1+reduction}(end+1,:) = [ytrue.VVrank{idx}(idx2),ytrue.VV{idx}(idx2)];
    end
end
figure;
hold on
for idx = 1: 4
    scatter(1 -todo{idx}(:,2), todo{idx}(:,1))
end
xlabel('1 - prob')
ylabel('rank')
set(gca, 'XScale','log', 'YScale', 'log');
title('agnostic')
legend('0', '-1', '-2', '-3');