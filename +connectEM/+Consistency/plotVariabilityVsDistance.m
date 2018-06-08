function plotVariabilityVsDistance( ...
        synT, synToSyn, pairConfig, varargin)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    opt = struct;
    opt.maxDistUm = [];
    opt.minDistUm = [];
    opt.info = Util.runInfo();
    opt = Util.modifyStruct(opt, varargin{:});
    
    curT = struct2table(rmfield(pairConfig, 'title'));
    curT.preAggloId = synT.preAggloId(curT.synIdPairs(:, 1));
    curT.postAggloId = synT.postAggloId(curT.synIdPairs(:, 1));

    curT.preDist = zeros(size(curT.preAggloId));
    curT.postDist = zeros(size(curT.postAggloId));
    for curIdx = 1:numel(curT.preDist)
        curPreAggloId = curT.preAggloId(curIdx);
        curPostAggloId = curT.postAggloId(curIdx);
        curSynIds = synT.id(curT.synIdPairs(curIdx, :));

        % Distance along axonal side
        curPreSynIds = synToSyn.axonSynIds{curPreAggloId};
       [~, curPreSynIds] = ismember(curSynIds, curPreSynIds);
        curPreDist = synToSyn.axonSynToSynDists{curPreAggloId};
        curPreDist = curPreDist(curPreSynIds(1), curPreSynIds(2));

        % Distance along axonal side
        curPostSynIds = synToSyn.dendSynIds{curPostAggloId};
       [~, curPostSynIds] = ismember(curSynIds, curPostSynIds);
        curPostDist = synToSyn.dendSynToSynDists{curPostAggloId};
        curPostDist = curPostDist(curPostSynIds(1), curPostSynIds(2));

        curT.preDist(curIdx) = curPreDist;
        curT.postDist(curIdx) = curPostDist;
    end

    curT.cv = synT.area(curT.synIdPairs);
    curT.cv = std(curT.cv, 0, 2) ./ mean(curT.cv, 2);
    
    if ~isempty(opt.maxDistUm)
        curPreT = curT(curT.preDist / 1E3 < opt.maxDistUm, :);
        curPostT = curT(curT.postDist / 1E3 < opt.maxDistUm, :);
    else
        curPreT = curT;
        curPostT = curT;
    end
    
    curFitPre = fit(curPreT.preDist / 1E3, curPreT.cv, 'poly1');
    curFitPost = fit(curPostT.postDist / 1E3, curPostT.cv, 'poly1');

    % Plot
    fig = figure();
    fig.Color = 'white';
    fig.Position(3:4) = [720, 480];

    % Presynaptic side
    ax = subplot(1, 2, 1);
    axis(ax, 'square');
    hold(ax, 'on');

    plot(ax, curPreT.preDist / 1E3, curPreT.cv, '.');
    plot(ax, ax.XLim, curFitPre(ax.XLim), 'Color', 'black');
    title(ax, {'Presynaptic', sprintf( ...
        'CV = %.2f + %.4fd', curFitPre.p2, curFitPre.p1)}, ...
        'FontWeight', 'normal', 'FontSize', 10);
    
    fitlm(curPreT.preDist / 1E3, curPreT.cv)

    % Postsynaptic side
    ax = subplot(1, 2, 2);
    axis(ax, 'square');
    hold(ax, 'on');

    plot(ax, curPostT.postDist / 1E3, curPostT.cv, '.');
    plot(ax, ax.XLim, curFitPost(ax.XLim), 'Color', 'black');
    title(ax, {'Postsynaptic', sprintf( ...
        'CV = %.2f + %.4fd', curFitPost.p2, curFitPost.p1)}, ...
        'FontWeight', 'normal', 'FontSize', 10);
    
    fitlm(curPostT.postDist / 1E3, curPostT.cv)

    axes = fig.Children;
    maxX = max(cat(2, axes.XLim));
    
    set(axes, ...
        'TickDir', 'out', ...
        'XLim', [0, maxX], ...
        'YLim', [0, sqrt(2)]);

    ax = axes(end);
    xlabel(ax, 'Intersynapse distance (µm)');
    ylabel(ax, 'Coefficient of variation');
    
    annotation(fig, ...
        'textbox', [0, 0.9, 1, 0.1], ...
        'String', { ...
            opt.info.filename; ...
            opt.info.git_repos{1}.hash; ...
            pairConfig.title}, ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'center');
end
