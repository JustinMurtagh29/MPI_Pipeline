function visualizeTasks(settings, missions, seg)

for task=1:length(missions)
    figure('Visible', 'off', 'Renderer', 'OpenGL');
    % Plot starting position & direction
    localPos = missions(task).start.position - settings.bbox(:,1)';
    plot3(localPos(2), localPos(1), localPos(3), 'xr', 'LineWidth', 5);
    hold on;
    direct = missions(task).start.direction;
    quiver3(localPos(2), localPos(1), localPos(3), direct(2), direct(1), direct(3), 0, 'r', 'LineWidth', 5);
    % Plot error center
    localCenter = missions(task).errorCenter - settings.bbox(:,1)';
    plot3(localCenter(2), localCenter(1), localCenter(3), 'og', 'LineWidth', 5);
    % Plot start object solid
    obj = seg == missions(task).start.id;
    issf = isosurface(obj, .1);
    k = patch(issf);
    set(k, 'FaceColor', 'g', 'EdgeColor', 'none', 'FaceAlpha', .8);
    % Plot other objects transparent
    colors = jet(length(missions(task).possibleEnds));
    for i=1:length(missions(task).possibleEnds)    
        obj = seg == missions(task).possibleEnds(i).id;
        issf = isosurface(obj, .1);
        l = patch(issf);
        set(l, 'FaceColor', colors(i,:), 'EdgeColor', 'none', 'FaceAlpha', .1);
    end
    view(3);
    daspect([28 28 11]);
    title(['Task: ' num2str(task)]);
    grid on;
    camlight('headlight');
    lighting phong;
    legend({'start position', 'start direction', 'error center', 'start object'});
    set(gcf, 'PaperPositionMode', 'manual', 'PaperType', 'A4', 'PaperUnits', 'normalized', 'PaperPosition', [0 0 1 1]);
    saveas(gcf, ['/zdata/manuel/sync/activeTraining/tasks/' num2str(missions(task).uniqueID, '%.3i') '.pdf']);
    close all;
end

end
