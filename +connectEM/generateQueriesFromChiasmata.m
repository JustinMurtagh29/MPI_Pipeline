function q = generateQueriesFromChiasmata(p, segmentMeta, input, outputFolder, options)
%connectEM.generateQueriesFromChiasmata(p,segmentMeta,output,'/gaba/scratch/kboerg/queriesMHchiasma/',struct('queryBoundingBoxSize',4000,'writeTasksToFile',1))
    q.pos = {};
    q.dir = {};
    q.order = {};
    voxelSize = p.raw.voxelSize;
    % Decide which positions to querry and calculate some statistics
    for idx = 1 : 10
        for idx2 = 1 : length(input{idx})
            if isfield(input{idx}{idx2}, 'position')
            for idx3 = 1 : length(input{idx}{idx2}.position)
                q.pos{end +1} = input{idx}{idx2}.position{idx3};
                q.dir{end +1} = -input{idx}{idx2}.direction{idx3}; %chiasma queries are reversed, see line dir{i} = bsxfun(@minus, pos{i}, nodesV(cc{i}(centerOfCC(i)),:));
                q.order{end + 1} = repmat(idx2 * 100 + idx, size(input{idx}{idx2}.position{idx3},1), 1);
            end
        end
        end
    end
    % Sort out all queries based on some heuristics
    % borderSize = round(options.datasetBorderExclusionSize ./ voxelSize);
    % bbox(:,1) = p.bbox(:,1) + borderSize';
    % bbox(:,2) = p.bbox(:,2) - borderSize';
    % q.outsideBBox = cellfun(@(x)any(bsxfun(@le, x, bbox(:,1)'),2) | any(bsxfun(@ge, x, bbox(:,2)'),2), q.pos, 'uni', 0);
    % % MH does not want these exclusion criteria anymore (for now)
    %q.tooSmall = cellfun(@(x)x < 5000, q.voxelCount);
    %q.tooShort = cellfun(@(x)x < 1500, q.lengthAlongPC1);
    %q.tooUnstraight = cellfun(@(x)x(1) < .7, q.varianceExplained);
    %q.exclude = arrayfun(@(x,y,z,t,u) x{:} | y | z | t, q.outsideBBox, q.tooShort, q.tooSmall, q.tooUnstraight, 'uni', 0);
    %q.exclude = arrayfun(@(x)x{:}, q.outsideBBox, 'uni', 0);
    % 'Write' problems to flight mode webKNOSSOS
    if options.writeTasksToFile
        extend = round(options.queryBoundingBoxSize ./ voxelSize);
        fid = fopen([outputFolder datestr(clock,30) '_flightTasks.txt'], 'w');
        for i=1:length(q.pos)
            for j=1:size(q.pos{i},1)
                [phi, theta, psi] = calculateEulerAngles(q.dir{i}(j,:), voxelSize);
                q.angles{i}(j,:) = [phi theta psi];
                minPos = q.pos{i}(j,:) - extend;
                sizeBbox = 2*extend;
                taskString = ['2012-09-28_ex145_07x2_ROI2017,56d6a7c6140000d81030701e,queriesMHchiasma,1,' ...
                    num2str(q.pos{i}(j,1)-1) ',' num2str(q.pos{i}(j,2)-1) ',' num2str(q.pos{i}(j,3)-1) ',' ...
                    num2str(phi) ',' num2str(theta) ',' num2str(psi) ',1,Connectomics department,' ...
                    num2str(minPos(1)) ',' num2str(minPos(2)) ',' num2str(minPos(3)) ',' ...
                    num2str(sizeBbox(1)) ',' num2str(sizeBbox(2)) ',' num2str(sizeBbox(3)) ',' 'queriesMHchiasma'];
                %if ~q.exclude{i}(j)
                    fprintf(fid, '%s\n', taskString);
                %send
            end
        end
        fclose(fid);
    end
end

function [p, d, vE, s, l] = determineQueryLocation(segmentMeta, agglo, voxelSize)
    % Determine query position, direction and edge based on PCA

    thesePoints = segmentMeta.point(:,agglo)';
    [coeff,score,latent] = pca(bsxfun(@times, thesePoints, voxelSize));
    [~,minIdx] = min(score(:,1));
    [~,maxIdx] = max(score(:,1));
    p(1,:) = thesePoints(minIdx,:);
    p(2,:) = thesePoints(maxIdx,:);
    % Transform direction (1.PC) back into voxel space for wK and normalize again
    d(1,:) = -coeff(:,1) .* 1./voxelSize';
    d(2,:) = coeff(:,1) .* 1./voxelSize';
    d(1,:) = d(1,:) ./ norm(d(1,:));
    d(2,:) = d(2,:) ./ norm(d(2,:));
    vE = latent./sum(latent);
    % Determine number of voxel for current agglo
    s = sum(segmentMeta.voxelCount(agglo));
    % Determine length of agglomeration along PC1
    l = sqrt(sum(((p(2,:) - p(1,:)).*voxelSize).^2));

end

function [phi, thetha, psi] = calculateEulerAngles(di, voxelSize)
    % Calculate angles (as deinfed in wK) in degrees from direction vector

    % Make sure direction is normalized
    di = di .* voxelSize;
    di = di ./ norm(di);
    % Calculate euler angles, from wk-paper -> RESCOPaaS
    eulerAngles = diffToEulerAngle(di);
    phi = round(eulerAngles(1));
    thetha = round(eulerAngles(2));
    psi = round(eulerAngles(3));

end
