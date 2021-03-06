function debugNewQueries(segmentMeta, agglos, q, outputFolder)
    % Visualize a set of agglomerates and queries generated by/from them

    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end

    theseNodes = {};
    theseNames = {};
    theseComments = {};
    for i=1:length(agglos)
        theseNodes{end+1} = segmentMeta.point(:,agglos{i})';
        theseNames{end+1} = [num2str(i, '%.3i') '_agglomerate'];
        theseComments{end+1} = {};
        for j=1
            % Write queries indicating direction (and comment to paste behing link for flight location)
            theseNodes{end+1}(1,:) = round(q.pos{i}(j,:));
            theseNodes{end}(2,:) = round(q.pos{i}(j,:) + 10 * q.dir{i}(j,:));
            theseNames{end+1} = [num2str(i, '%.3i') '_query ' num2str(j)];
            %if q.exclude{i}(j)
            %    theseNames{end} = [theseNames{end} '_excluded'];
            %end
            theseComments{end+1}{1} = ['Pos: ' num2str(q.pos{i}(j,1)) ',' num2str(q.pos{i}(j,2)) ',' num2str(q.pos{i}(j,3)) ',' ...
                ' Angles: ' num2str(q.angles{i}(j,1)) ',' num2str(q.angles{i}(j,2)) ',' num2str(q.angles{i}(j,3))];
        end
    end
    connectEM.generateSkeletonFromNodes( ...
        fullfile(outputFolder,[datestr(clock, 30) '_skel.nml']), ...
        theseNodes, theseNames, theseComments);
end

