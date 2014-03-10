function visualizeTasksNew(seg, v, parameter, kdbParameter, id, i, colors)
seg = permute(seg,[2 1 3]);

for vp =1:length(v)
    possibleEnds = [v(vp).missions(i).possibleEnds(:).id];
    idRenderedEnd = find(possibleEnds == v(vp).missions(i).end.id);
    resortedIdx = 1:length(v(vp).missions(i).possibleEnds);
    resortedIdx(1) = idRenderedEnd(1);
    resortedIdx(idRenderedEnd(1)) = 1;
    v(vp).missions(i).possibleEnds = v(vp).missions(i).possibleEnds(resortedIdx);
    % Rendering Length
    startFF(i) = v(vp).missions(i).start.firstFrame;
    startLF(i) = v(vp).missions(i).start.lastFrame;
    for j=1:length(v(vp).missions(i).possibleEnds)
        endFF{i,j} = v(vp).missions(i).possibleEnds(j).firstFrame;
        endLF{i,j} = v(vp).missions(i).possibleEnds(j).lastFrame;
    end
    figure('Units', 'centimeters', 'Position', [0 0 29.7 21], ...
    'PaperType', 'A4', 'PaperUnits', 'centimeters', 'PaperPosition', [0 0 29.7 21], 'PaperOrientation', 'portrait');
    % Plot CoM starting object
    p = v(vp).missions(i).start.CoM .* kdbParameter.settings.scale;
    h1 = plot3(p(1), p(2), p(3), 'xr', 'MarkerSize', 10, 'LineWidth', 3);
    hold on;
    % Plot CoM end object
    p = v(vp).missions(i).end.CoM .* kdbParameter.settings.scale;
    h2 = plot3(p(1), p(2), p(3), 'xg', 'MarkerSize', 10, 'LineWidth', 3);
    % Plot error center
    p = v(vp).missions(i).errorCenter .* kdbParameter.settings.scale;
    h3 = plot3(p(1), p(2), p(3), 'xy', 'MarkerSize', 10, 'LineWidth', 3);
    % Plot forward arrow according to rendering length
    d = (v(vp).missions(i).start.lastFrame - v(vp).missions(i).start.firstFrame) .* v(vp).missions(i).start.direction;
    p = v(vp).missions(i).errorCenter .* kdbParameter.settings.scale + v(vp).missions(i).start.firstFrame*v(vp).missions(i).start.direction;  
    h4 = quiver3(p(1),p(2),p(3),d(1),d(2),d(3), 0, 'r', 'LineWidth', 3);
    % Plot backward arrow according to rendering length
    d = (v(vp).missions(i).possibleEnds(1).lastFrame - v(vp).missions(i).possibleEnds(1).firstFrame) .* v(vp).missions(i).start.direction;
    p = v(vp).missions(i).errorCenter .* kdbParameter.settings.scale + v(vp).missions(i).possibleEnds(1).firstFrame*v(vp).missions(i).start.direction; 
    h5 = quiver3(p(1),p(2),p(3),d(1),d(2),d(3), 0, 'g', 'LineWidth', 3);
    % Plot arrow according width of viewport
    d = 11.28 .* v(vp).width/2 .* v(vp).missions(i).start.orth1;
    p = v(vp).missions(i).errorCenter .* kdbParameter.settings.scale;  
    h6 = quiver3(p(1),p(2),p(3),d(1),d(2),d(3), 0, 'c', 'LineWidth', 3);
    % Plot arrow according height of viewport
    d = 11.28 .* v(vp).height/2 .* v(vp).missions(i).start.orth2;
    p = v(vp).missions(i).errorCenter .* kdbParameter.settings.scale;  
    quiver3(p(1),p(2),p(3),d(1),d(2),d(3), 0, 'c', 'LineWidth', 3);
%         % Plot inner bounding box for graph construction
%         bbox = parameter.bboxSmall' .* repmat(kdbParameter.settings.scale, 2, 1);
%         h7 = plot3([bbox(1,1) bbox(2,1)], [bbox(1,2) bbox(1,2)], [bbox(1,3) bbox(1,3)], ':m', 'LineWidth', 3);
%         plot3([bbox(1,1) bbox(2,1)], [bbox(2,2) bbox(2,2)], [bbox(1,3) bbox(1,3)], ':m', 'LineWidth', 3);
%         plot3([bbox(1,1) bbox(2,1)], [bbox(1,2) bbox(1,2)], [bbox(2,3) bbox(2,3)], ':m', 'LineWidth', 3);
%         plot3([bbox(1,1) bbox(2,1)], [bbox(2,2) bbox(2,2)], [bbox(2,3) bbox(2,3)], ':m', 'LineWidth', 3);
%         plot3([bbox(2,1) bbox(2,1)], [bbox(1,2) bbox(2,2)], [bbox(1,3) bbox(1,3)], ':m', 'LineWidth', 3);
%         plot3([bbox(1,1) bbox(1,1)], [bbox(1,2) bbox(2,2)], [bbox(2,3) bbox(2,3)], ':m', 'LineWidth', 3);
%         plot3([bbox(2,1) bbox(2,1)], [bbox(1,2) bbox(2,2)], [bbox(2,3) bbox(2,3)], ':m', 'LineWidth', 3);
%         plot3([bbox(1,1) bbox(1,1)], [bbox(1,2) bbox(2,2)], [bbox(1,3) bbox(1,3)], ':m', 'LineWidth', 3);
%         plot3([bbox(1,1) bbox(1,1)], [bbox(1,2) bbox(1,2)], [bbox(1,3) bbox(2,3)], ':m', 'LineWidth', 3);
%         plot3([bbox(1,1) bbox(1,1)], [bbox(2,2) bbox(2,2)], [bbox(1,3) bbox(2,3)], ':m', 'LineWidth', 3);
%         plot3([bbox(2,1) bbox(2,1)], [bbox(1,2) bbox(1,2)], [bbox(1,3) bbox(2,3)], ':m', 'LineWidth', 3);
%         plot3([bbox(2,1) bbox(2,1)], [bbox(2,2) bbox(2,2)], [bbox(1,3) bbox(2,3)], ':m', 'LineWidth', 3);
    % Plot start object rather solid
    obj = seg == v(vp).missions(i).start.id;
    issf = isosurface(obj, .1);
    issf.vertices = bsxfun(@times, bsxfun(@plus, issf.vertices, parameter.bboxBig(:,1)'), kdbParameter.settings.scale);
    h8 = patch(issf);
    set(h8, 'FaceColor', 'r', 'EdgeColor', 'none', 'FaceAlpha', .5);
    % Plot end object rather solid
    obj = seg == v(vp).missions(i).end.id;
    issf = isosurface(obj, .1);
    issf.vertices = bsxfun(@times, bsxfun(@plus, issf.vertices, parameter.bboxBig(:,1)'), kdbParameter.settings.scale);
    h9 = patch(issf);
    set(h9, 'FaceColor', 'g', 'EdgeColor', 'none', 'FaceAlpha', .5);
    % Plot other objects
%         for j=2:length(v(vp).missions(i).possibleEnds)
%                 obj = seg == v(vp).missions(i).possibleEnds(j).id;
%                 issf = isosurface(obj, .1);
%                 issf.vertices = bsxfun(@times, bsxfun(@plus, issf.vertices, parameter.bboxBig(:,1)'), kdbParameter.settings.scale);
%                 l(j) = patch(issf);
%                 set(l(j), 'FaceColor', colors(1+j,:), 'EdgeColor', 'none', 'FaceAlpha', .5);
%         end
    % Plot bbox used for local computation of direction/ rendering length
    bbox = [(v(vp).missions(i).errorCenter - kdbParameter.border).* kdbParameter.settings.scale; ...
        (v(vp).missions(i).errorCenter + kdbParameter.border).* kdbParameter.settings.scale];
    h10 = plot3([bbox(1,1) bbox(2,1)], [bbox(1,2) bbox(1,2)], [bbox(1,3) bbox(1,3)], ':m', 'LineWidth', 3);
    plot3([bbox(1,1) bbox(2,1)], [bbox(2,2) bbox(2,2)], [bbox(1,3) bbox(1,3)], ':m', 'LineWidth', 3);
    plot3([bbox(1,1) bbox(2,1)], [bbox(1,2) bbox(1,2)], [bbox(2,3) bbox(2,3)], ':m', 'LineWidth', 3);
    plot3([bbox(1,1) bbox(2,1)], [bbox(2,2) bbox(2,2)], [bbox(2,3) bbox(2,3)], ':m', 'LineWidth', 3);
    plot3([bbox(2,1) bbox(2,1)], [bbox(1,2) bbox(2,2)], [bbox(1,3) bbox(1,3)], ':m', 'LineWidth', 3);
    plot3([bbox(1,1) bbox(1,1)], [bbox(1,2) bbox(2,2)], [bbox(2,3) bbox(2,3)], ':m', 'LineWidth', 3);
    plot3([bbox(2,1) bbox(2,1)], [bbox(1,2) bbox(2,2)], [bbox(2,3) bbox(2,3)], ':m', 'LineWidth', 3);
    plot3([bbox(1,1) bbox(1,1)], [bbox(1,2) bbox(2,2)], [bbox(1,3) bbox(1,3)], ':m', 'LineWidth', 3);
    plot3([bbox(1,1) bbox(1,1)], [bbox(1,2) bbox(1,2)], [bbox(1,3) bbox(2,3)], ':m', 'LineWidth', 3);
    plot3([bbox(1,1) bbox(1,1)], [bbox(2,2) bbox(2,2)], [bbox(1,3) bbox(2,3)], ':m', 'LineWidth', 3);
    plot3([bbox(2,1) bbox(2,1)], [bbox(1,2) bbox(1,2)], [bbox(1,3) bbox(2,3)], ':m', 'LineWidth', 3);
    plot3([bbox(2,1) bbox(2,1)], [bbox(2,2) bbox(2,2)], [bbox(1,3) bbox(2,3)], ':m', 'LineWidth', 3);
    % Plot viewport
    step = [-1 1];
    for x=1:2
        for y=1:2
            for z=1:2
                edgePos{x, y, z} = v(vp).missions(i).errorCenter .* kdbParameter.settings.scale + ...
                    step(x) * 11.28 .* v(vp).width/2 .* v(vp).missions(i).start.orth1 + ...
                    step(y) * 11.28 .* v(vp).height/2 .* v(vp).missions(i).start.orth2 + ...
                    step(z) * 11.28 .* 25 .* v(vp).missions(i).start.direction;
            end
        end
    end
    h11 = plot3([edgePos{1,1,1}(1) edgePos{1,1,2}(1)], [edgePos{1,1,1}(2) edgePos{1,1,2}(2)], [edgePos{1,1,1}(3) edgePos{1,1,2}(3)], ':c', 'LineWidth', 3);
    plot3([edgePos{1,1,1}(1) edgePos{1,2,1}(1)], [edgePos{1,1,1}(2) edgePos{1,2,1}(2)], [edgePos{1,1,1}(3) edgePos{1,2,1}(3)], ':c', 'LineWidth', 3);
    plot3([edgePos{1,1,1}(1) edgePos{2,1,1}(1)], [edgePos{1,1,1}(2) edgePos{2,1,1}(2)], [edgePos{1,1,1}(3) edgePos{2,1,1}(3)], ':c', 'LineWidth', 3);    
    plot3([edgePos{2,2,2}(1) edgePos{1,2,2}(1)], [edgePos{2,2,2}(2) edgePos{1,2,2}(2)], [edgePos{2,2,2}(3) edgePos{1,2,2}(3)], ':c', 'LineWidth', 3);
    plot3([edgePos{2,2,2}(1) edgePos{2,2,1}(1)], [edgePos{2,2,2}(2) edgePos{2,2,1}(2)], [edgePos{2,2,2}(3) edgePos{2,2,1}(3)], ':c', 'LineWidth', 3);
    plot3([edgePos{2,2,2}(1) edgePos{2,1,2}(1)], [edgePos{2,2,2}(2) edgePos{2,1,2}(2)], [edgePos{2,2,2}(3) edgePos{2,1,2}(3)], ':c', 'LineWidth', 3);
    plot3([edgePos{2,2,1}(1) edgePos{2,1,1}(1)], [edgePos{2,2,1}(2) edgePos{2,1,1}(2)], [edgePos{2,2,1}(3) edgePos{2,1,1}(3)], ':c', 'LineWidth', 3);
    plot3([edgePos{2,2,1}(1) edgePos{1,2,1}(1)], [edgePos{2,2,1}(2) edgePos{1,2,1}(2)], [edgePos{2,2,1}(3) edgePos{1,2,1}(3)], ':c', 'LineWidth', 3);
    plot3([edgePos{1,2,2}(1) edgePos{1,1,2}(1)], [edgePos{1,2,2}(2) edgePos{1,1,2}(2)], [edgePos{1,2,2}(3) edgePos{1,1,2}(3)], ':c', 'LineWidth', 3);
    plot3([edgePos{1,2,2}(1) edgePos{1,2,1}(1)], [edgePos{1,2,2}(2) edgePos{1,2,1}(2)], [edgePos{1,2,2}(3) edgePos{1,2,1}(3)], ':c', 'LineWidth', 3);
    plot3([edgePos{2,1,2}(1) edgePos{2,1,1}(1)], [edgePos{2,1,2}(2) edgePos{2,1,1}(2)], [edgePos{2,1,2}(3) edgePos{2,1,1}(3)], ':c', 'LineWidth', 3);
    plot3([edgePos{2,1,2}(1) edgePos{1,1,2}(1)], [edgePos{2,1,2}(2) edgePos{1,1,2}(2)], [edgePos{2,1,2}(3) edgePos{1,1,2}(3)], ':c', 'LineWidth', 3);
    view(3);
    title(['Mission ID: ' num2str(i) ', Game: ' v(vp).name ', Error center: ' num2str(v(vp).missions(i).errorCenter)]);
    grid on;
    camlight('headlight');
    lighting phong;
    legend([h1 h2 h3 h4 h5 h6 h8 h9 h10 h11], 'CoM start object', 'CoM end object', 'CoM border', 'rendering length start', 'rendering length end', 'ONB', 'start object', 'end object', 'local bbox', 'viewport');
    axis equal;
    saveas(gcf, ['C:\Users\mberning\Desktop\problemInspector\' id v(vp).name num2str(i, '%.3i') 'isosurface.fig']);
    set(gcf,'PaperSize',fliplr(get(gcf,'PaperSize')));
    print(gcf, ['C:\Users\mberning\Desktop\problemInspector\' id v(vp).name num2str(i, '%.3i') 'isosurface.pdf'], '-dpdf', '-r300');
    close all;

    %% Second: rendering length & probability visualization
    figure('Units', 'centimeters', 'Position', [0 0 29.7 21], ...
    'PaperType', 'A4', 'PaperUnits', 'centimeters', 'PaperPosition', [0 0 29.7 21], 'PaperOrientation', 'landscape');
    % rendering length
    subplot(1,2,1);
    hold on;
    plot([0 0], [startFF(i) startLF(i)], 'Color', colors(1,:), 'LineWidth', 3);
    labels{1} = num2str(v(vp).missions(i).start.id);
    for j=1:length(v(vp).missions(i).possibleEnds)
        if ~isempty(endFF{i,j})
            plot([j j], [endFF{i,j} endLF{i,j}], 'Color', colors(1+j,:), 'LineWidth', 3);
            labels{1+j} = num2str(v(vp).missions(i).possibleEnds(j).id);
        end
        plot(j, v(vp).missions(i).possibleEnds(j).touchFrame, 'x', 'Color', colors(1+j,:));
    end
    title(['Task ' num2str(i) ': original, rendering length']);
    xlabel('Object IDs');
    xlim([-0.5 length(v(vp).missions(i).possibleEnds)+0.5]);
    set(gca, 'Xtick', 0:1:length(v(vp).missions(i).possibleEnds));
    set(gca, 'XTickLabel', labels);
    ylabel('distance along problem direction [nm]');
    % probability
    subplot(1,2,2);
    hold on;
    x = [v(vp).missions(i).possibleEnds(:).probability]' * 100;
    h = bar(x);
    ch = get(h,'Children'); %get children of the bar group
    fvd = get(ch,'Faces'); %get faces data
    fvcd = get(ch,'FaceVertexCData'); %get face vertex cdata
    for k = 1:length(x)
        row = k;
        fvcd(fvd(row,:)) = k; %adjust the face vertex cdata to be that of the row
    end
    set(ch,'FaceVertexCData',fvcd) %set to new face vertex cdata
    colormap(colors(2:length(v(vp).missions(i).possibleEnds),:));
    set(gca, 'Clim', [1 length(v(vp).missions(i).possibleEnds)]);
    labels = labels(2:end);
    title(['Task ' num2str(i) ': original, probabilities']);
    set(gca, 'Xtick', 1:length(v(vp).missions(i).possibleEnds));
    set(gca, 'XTickLabel', labels);
    xlim([0.5 length(v(vp).missions(i).possibleEnds)+0.5]);
    xlabel('Object IDs');
    ylabel('probability [%]');    
    saveas(gcf, ['C:\Users\mberning\Desktop\problemInspector\' id v(vp).name num2str(i, '%.3i') 'statistics.pdf']);
    close all;

end

