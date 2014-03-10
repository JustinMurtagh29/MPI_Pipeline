function visualizeEdgeFeaturesForErrors( weights, label, errorIdx, merge, discard, isConnected)
% Visualizierungs-Settings
% colors = {'g' 'r' 'y' 'g' 'r' 'y' 'g' 'r' 'y'};
% lineStyle = {':' ':' ':' '-' '-' '-' '--' '--' '--'};
colors = {'g' 'r' 'b' 'y'};
lineStyle = {'-' '-' '-' '-'};
% Training Daten aufsplitten
if merge(errorIdx)
    tp = weights(merge & (isConnected == 1),:);
    fp = weights(merge & (isConnected == -1),:);
    tn = weights(~merge & (isConnected == -1),:);
    fn = weights(~merge & (isConnected == 1),:);
end
if discard(errorIdx)
    tp = weights(discard & (isConnected == -1),:);
    fp = weights(discard & (isConnected == 1),:);
    tn = weights(~discard & (isConnected == 1),:);
    fn = weights(~discard & (isConnected == -1),:);
end
% tm = weights(merge & (isConnected == 1),:);
% fm = weights(merge & (isConnected == -1),:);
% um = weights(merge & (isConnected == 0),:);
% td = weights(discard & (isConnected == -1),:);
% fd = weights(discard & (isConnected == 1),:);
% ud = weights(discard & (isConnected == 0),:);
% pq = weights(~discard & ~merge & (isConnected == 1),:);
% nq = weights(~discard & ~merge & (isConnected == -1),:);
% uq = weights(~discard & ~merge & (isConnected == 0),:);
for idx = 1:215
    ax = subplot(22,10,idx);
    % Bin Vektor berechnen
%     limits = [min([tm(:,idx); fm(:,idx); um(:,idx); td(:,idx); fd(:,idx); ud(:,idx); pq(:,idx); nq(:,idx); uq(:,idx);]) ...
%         max([tm(:,idx); fm(:,idx); um(:,idx); td(:,idx); fd(:,idx); ud(:,idx); pq(:,idx); nq(:,idx); uq(:,idx);])];
    limits = [min([tp(:,idx); fp(:,idx); tn(:,idx); fn(:,idx)]) ...
        max([tp(:,idx); fp(:,idx); tn(:,idx); fn(:,idx)])];    
    binSize = (limits(2) - limits(1))/20;
    x = limits(1)+binSize/2:binSize:limits(2)-binSize/2;
    % True Positives
    a = hist(tp(:,idx),x);
    a = a ./ sum(a);
    % False Positives
    b = hist(fp(:,idx),x);
    b = b ./ sum(b);
    % True Negatives
    c = hist(tn(:,idx),x);
    c = c ./ sum(c);
    % False Negatives
    d = hist(fn(:,idx),x);
    d = d ./ sum(d);
%     % True Merger
%     a = hist(tm(:,idx),x);
%     a = a ./ sum(a);
%     % False Merger
%     b = hist(fm(:,idx),x);
%     b = b ./ sum(b);
%     % Unlabeled Merger
%     c = hist(um(:,idx),x);
%     c = c ./ sum(c);
%     % True Discards
%     d = hist(td(:,idx),x);
%     d = d ./ sum(d);
%     % False Discards
%     e = hist(td(:,idx),x);
%     e = e ./ sum(e);
%     % Unlabeled Discards
%     f = hist(ud(:,idx),x);
%     f = f ./ sum(f);
%     % Positive Querries
%     g = hist(pq(:,idx),x);
%     g = g ./ sum(g);
%     % Negative Querries
%     h = hist(nq(:,idx),x);
%     h = h ./ sum(h);
%     % Unlabeled Querries
%     i = hist(uq(:,idx),x);
%     i = i ./ sum(i);
    % Stairs-Plot
    fh = stairs(x, [a; b; c; d]');
%     fh = stairs(x, [a; b; c; d; e; f; g; h; i]');
    for j=1:length(fh)
        set(fh(j),'Color',colors{j}, 'LineStyle', lineStyle{j}, 'LineWidth', .5);
    end
    hold on;
    ylims = get(gca, 'YLim');
    plot([weights(errorIdx,idx) weights(errorIdx,idx)], ylims, 'k', 'LineWidth', .5);
    xlim(limits);
    %text(0.5, 0.5, sprintf(strrep(label{idx}, ' ', '\n')));
    set(ax, 'visible', 'off');
    pos = get(ax, 'Position');
    pos(1:2) = pos(1:2) - pos(3:4)/7;
    pos(3:4) = pos(3:4) + 2*pos(3:4)/7;
    set(ax, 'Position', pos);
    %legend('true merger', 'false merger', 'unlabeled merger', 'true discard', 'false discard', 'unlabeled discard', 'positive query', 'negative query', 'unlabeled query', 'value for error');
    %legend('true positive', 'false positive', 'true negative', 'false negative', 'value for error');
end

end
