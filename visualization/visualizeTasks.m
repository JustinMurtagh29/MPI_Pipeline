function visualizeTasks(seg, edges, p, edgesOld, pOld, id, problemNumbers)

edgesToShow = edges(problemNumbers,:);

for i=1:size(edgesToShow,1)
    idx1 = edgesToShow(i,1);
    idx2 = edgesToShow(i,2);
    obj1 = seg == idx1;
    obj2 = seg == idx2;
    % Determine CoM first object
    obj1Prop = regionprops(obj1, 'Centroid', 'Area');
    [~, idxx] = max([obj1Prop(:).Area]);
    obj1Prop = obj1Prop(idxx);
    % Determine CoM second object
    obj2Prop = regionprops(obj2, 'Centroid', 'Area');
    [~, idxx] = max([obj2Prop(:).Area]);
    obj2Prop = obj2Prop(idxx);
    % Determine CoM border between objects (largest area of touch)
    border = bwdist(obj1) <= 1 & bwdist(obj2) <= 1;
    borderProp = regionprops(border, 'Centroid', 'Area');
    [~, idx] = max([borderProp(:).Area]);
    % Mission definitions
    mission.errorCenter = round(borderProp(idx).Centroid([2 1 3]));
    mission.start.direction = obj2Prop.Centroid([2 1 3]) - obj1Prop.Centroid([2 1 3]);
    mission.start.direction = mission.start.direction / norm(mission.start.direction);
    mission.start.position = round(mission.errorCenter - 15 .* [1 1 0.5] .*  mission.start.direction);
    mission.start.id = idx1;
    [row, col] = find(edgesOld == idx1);
    col(col == 1) = 3;
    col = col - 1;
    for k=1:length(row)
        mission.possibleEnds(k).id = edgesOld(row(k),col(k));
        mission.possibleEnds(k).probability = pOld(row(k));
    end
    mission.difficulty = 1 - max([mission.possibleEnds(:).probability]);
    %% TILL HERE WRITEKNOWLEDGEDB.m, now visualization
    % Plot starting position
    localPos = obj1Prop.Centroid([2 1 3]);
    plot3(localPos(2), localPos(1), localPos(3), 'ob', 'MarkerSize', 7, 'LineWidth', 7);
    hold on;
    localPos = obj2Prop.Centroid([2 1 3]);
    plot3(localPos(2), localPos(1), localPos(3), 'og', 'MarkerSize', 7, 'LineWidth', 7);
    % Plot error center & direction
    localCenter = mission.errorCenter;
    plot3(localCenter(2), localCenter(1), localCenter(3), 'or','MarkerSize', 7, 'LineWidth', 7);
    direct = mission.start.direction;
    quiver3(localCenter(2), localCenter(1), localCenter(3), 100*direct(2), 100*direct(1), 100*direct(3), 0, 'g', 'LineWidth', 5);
    quiver3(localCenter(2), localCenter(1), localCenter(3), -100*direct(2), -100*direct(1), -100*direct(3), 0, 'b', 'LineWidth', 5);
    % Plot start object solid
    obj = seg == mission.start.id;
    issf = isosurface(obj, .1);
    k = patch(issf);
    set(k, 'FaceColor', 'b', 'EdgeColor', 'none', 'FaceAlpha', .1);
    % Plot other objects transparent
    for j=1:length(mission.possibleEnds)    
        obj = seg == mission.possibleEnds(j).id;
        issf = isosurface(obj, .1);
        l = patch(issf);
        color = [1-mission.possibleEnds(j).probability mission.possibleEnds(j).probability 0];
        set(l, 'FaceColor', color, 'EdgeColor', 'none', 'FaceAlpha', .1);
    end
    view(3);
    daspect([28 28 11]);
    title(['Task ' num2str(problemNumbers(i)) ': SizeStart ' num2str(obj1Prop.Area) ', SizeEnd ' num2str(obj2Prop.Area) ', SizeBorder ' num2str(borderProp(idx).Area)]);
    grid on;
    camlight('headlight');
    lighting phong;
    daspect([28 28 11.28]);
    legend('CoM start object', 'CoM end object', 'CoM error center', 'forward direction', 'backward direction', 'start object');
    saveas(gcf, ['C:\sync\problemInspector\' id num2str(problemNumbers(i), '%.3i') '.fig']);
    close all;
    clear mission;
end

end

