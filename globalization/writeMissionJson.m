function writeMissionJson(p, problems, isControl, isTutorial)

    % Load global counter variables
    if exist(p.kdb.counter, 'file')
        load(p.kdb.counter);
    else
        counter.tasks = 0;
    end

    % put numbers in setParameterSettings (needs refactoring soon)
    kdbBorder = [225 225 90];
    missionsFolder = [p.saveFolder 'missions/'];
    if ~exist(missionsFolder, 'dir')
        mkdir(missionsFolder);
    end
    visualizationFolder = [missionsFolder 'visualization/'];
    if ~exist(visualizationFolder, 'dir')
        mkdir(visualizationFolder);
    end

    % Metadata fuer Tasks generieren
    % Iterations of seeded processes
    for it=1:length(problems)
        % Current node
        for no=1:length(problems(it).node)
            % All neighbours of current node
            for ne=1:length(problems(it).node(no).CoM)
                % all borders between current node and respective neighboue
                for bo=1:size(problems(it).node(no).CoM{ne},1)
                    tic;
                    % Start and End segmentation ID for this problem
                    startId = problems(it).node(no).self.id; 
                    endId = problems(it).node(no).neighbours(ne).id;
                    % cut out local stack for speed purposes (and for having local CoM/direction)
                    errorCenter = problems(it).node(no).CoM{ne}(bo,:);
                    bbox = bsxfun(@plus, errorCenter, [-kdbBorder; kdbBorder]);
                    segLocal = readKnossosRoi(p.seg.root, p.seg.prefix, bbox', 'uint32', '', 'raw');
                    % Calculate object properties of start and end segment within local stack
                    obj1Prop = calculateObjectProperties(segLocal, startId); 
                    obj2Prop = calculateObjectProperties(segLocal, endId); 
                    % Calculate mission direction (normalized) in nm space
                    direction = (obj2Prop.Centroid - obj1Prop.Centroid) .* p.kdb.settings.scale;
                    direction = direction ./ norm(direction);
                    % First: Find a orthogonal vector & 2 more vec for ONB (uniformly sampled from vectors in orthogonal plane)
                    [randomOrthVec, orth1, orth2] = findRandomOrthogonal(direction);
                    % Calculate extend of start object wrt to error center in nm 
                    pixelPosWrtErrorCenter = bsxfun(@times, bsxfun(@minus, obj1Prop.PixelList, kdbBorder), par.settings.scale);
                    % Determine pixel list all objects
                    objAllProps = regionprops(segLocal, 'PixelList');
                    % Find objects ID of objects neighbouring object 1
                    [row, col] = find(edgesOld == edges(i,1));
                    col(col == 1) = 3;
                    col = col - 1;
                    for j=1:length(v)
                        % Set mission id
                        v(j).missions(i).missionId = taskId + i;
                        % Set control problem flag
                        v(j).missions(i).isControl = isControl;
                        v(j).missions(i).isTutorial = isTutorial;
                        % Start and end segment used for rendering (section-local)
                        v(j).missions(i).start.id = startId;
                        v(j).missions(i).end.id = endId;
                        % error center in voxel (global)
                        v(j).missions(i).errorCenter = errorCenter;
                        % Calculate CoM of objects in voxel
                        v(j).missions(i).start.CoM = obj1Prop.Centroid + errorCenter - kdbBorder;
                        v(j).missions(i).end.CoM = obj2Prop.Centroid + errorCenter - kdbBorder;
                        % Rotated or original direction
                        if v(j).isRotated
                            v(j).missions(i).start.direction = rotate(direction, randomOrthVec, rand(1)*v(j).openingAngle);
                            [~, v(j).missions(i).start.orth1, v(j).missions(i).start.orth2] = findRandomOrthogonal(v(j).missions(i).start.direction);
                        else
                            v(j).missions(i).start.direction = direction;
                            v(j).missions(i).start.orth1 = orth1;
                            v(j).missions(i).start.orth2 = orth2;
                        end
                        % Calculate rendering length
                        [v(j).missions(i).start.firstFrame v(j).missions(i).start.lastFrame] = ...
                            calcRL(pixelPosWrtErrorCenter, v(j).missions(i).start.direction, v(j).missions(i).start.orth1, v(j).missions(i).start.orth2, ...
                            v(j).width/2*par.settings.scale(1), v(j).height/2*par.settings.scale(2));
                        % Iterate over possible ends
                        currentEnd = 1;
                        for k=1:length(row)
                            if ~isempty(find(segLocal == edgesOld(row(k),col(k))))
                                % Write id of connected object and probability of connection to struct
                                v(j).missions(i).possibleEnds(currentEnd).id = edgesOld(row(k),col(k));
                                v(j).missions(i).possibleEnds(currentEnd).probability = pOld(row(k));
                                % Calculate rendering length
                                temp =  bsxfun(@times, bsxfun(@minus, objAllProps(v(j).missions(i).possibleEnds(currentEnd).id).PixelList(:,[2 1 3]), par.kdbBorder), par.settings.scale);
                                [v(j).missions(i).possibleEnds(currentEnd).firstFrame v(j).missions(i).possibleEnds(currentEnd).lastFrame] = ...
                                    calcRL(temp, v(j).missions(i).start.direction, v(j).missions(i).start.orth1, v(j).missions(i).start.orth2, ...
                                    v(j).width/2*par.settings.scale(1), v(j).height/2*par.settings.scale(2));
                                temp =  (borderOld(row(k)).Centroid([2 1 3]) - errorCenter) .* par.settings.scale;
                                v(j).missions(i).possibleEnds(currentEnd).touchFrame = temp * v(j).missions(i).start.direction';
                                % Increase counter for possible ends
                                currentEnd = currentEnd + 1;
                            end
                        end
                    end
                end
            end
        end
        % Some bloated code just to have an updating progress withoit flooding command line
        t = toc;
        if i==1
            fprintf(['\nMission just written (total ' num2str(size(edges,1)) '): ']);
        else
            for ch=1:charLastIter
                fprintf('\b');
            end
        end
        charLastIter = fprintf([num2str(i, '%.4i') ', time last mission: ' num2str(t, '%4.2f') ' seconds']);
    end
    fprintf('\n');

    % Save global counter variables
    save(p.counter, 'counter');

    % Write mission json (named according to input file)
    missionStruct.sectionId = 'section1';
    missionStruct.viewports = v;
    writeJson([missionsFolder 'missions.json'], missionStruct);

    % save everything to synced folder for problem inspection if flag is set
    if par.saveForProblemInspector
        % call some function I guess 
    end

end

function obj1Prop = calculateObjectProperties(segLocal, segId)
    obj1 = segLocal == segId;
    obj1Prop = regionprops(obj1, 'Centroid', 'Area', 'PixelList');
    % Invert x-y coordinated of Centroid
    obj1Prop.Centroid = obj1Prop.Centroid([2 1 3]);
    obj1Prop.PixelList = obj1Prop.PixelList(:, [2 1 3]);
    % This is necessary if object is split by by cutting out local bbox
    if length(obj1Prop) > 1
        comWrtEc = sum(bsxfun(@times, bsxfun(@minus, cat(1,obj1Prop(:).Centroid), par.kdbBorder), par.settings.scale) .^ 2,2);
        [~, comIdx] = min(comWrtEc);
        obj1Prop = obj1Prop(comIdx);
    end
end

function [randomVec, orth1, orth2] = findRandomOrthogonal(v)
    % returns random orthogonal vector and and the 2nd and 3rd vector for ONB 
    v = v ./ norm(v);
    if all(v == [0 0 1])
        orth1 = [1 0 -1].*v;  
    else
        orth1 = [1 -1 0].*v;
    end
    orth1 = orth1 ./ norm(orth1);
    orth2 = cross(v, orth1); 
    mix = 2 * (rand(2,1) -0.5);
    randomVec = mix(1)*orth1+mix(2)*orth2;
    randomVec = randomVec ./ norm(randomVec);
end

function rotatedVec = rotate(vector,axis,phi)
    % vector 1 x 3, axis 1 x 3, not checked in here
    L0=axis/norm(axis);
    cphi=cosd(phi);
    rotatedVec=vector*cphi+(vector*L0')*(1-cphi)*L0+cross(L0,vector)*sind(phi);
end

function [ff, lf] = calcRL(pixelPosWrtErrorCenter, direction, orth1, orth2, width, height);
    objProjectionOntoDirection = pixelPosWrtErrorCenter * direction';
    objProjectionOntoOrth1 = pixelPosWrtErrorCenter * orth1';
    objProjectionOntoOrth2 = pixelPosWrtErrorCenter * orth2';
    outOfViewport = abs(objProjectionOntoOrth1) > width | abs(objProjectionOntoOrth2) > height;
    objProjectionOntoDirection = objProjectionOntoDirection(~outOfViewport);
    ff = min(objProjectionOntoDirection);
    lf = max(objProjectionOntoDirection);
end

