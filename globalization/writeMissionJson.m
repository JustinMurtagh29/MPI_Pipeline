function writeMissionJson(par, problems)

    % Load global counter variables
    if exist(par.counter, 'file')
        load(par.counter);
    else
        counter.taskId = 0;
    end

    % Metadata fuer Tasks generieren
    for it=1:length(problems)
        for no=1:length(problems(it).node)
        tic;
        % cut out local stack for speed purposes (and for having local CoM/direction)
        errorCenter = problems(it).CoM;
        bbox = errorCenter - [200 200 200]
        segLocal = seg(errorCenter(1)-par.kdbBorder(1):errorCenter(1)+par.kdbBorder(1),...
            errorCenter(2)-par.kdbBorder(2):errorCenter(2)+par.kdbBorder(2),...
            errorCenter(3)-par.kdbBorder(3):errorCenter(3)+par.kdbBorder(3));
        % Logical array of objects in small bounding box
        obj1 = segLocal == edges(i,1);
        obj2 = segLocal == edges(i,2);
        % Determine CoM of segments
        obj1Prop = regionprops(obj1, 'Centroid', 'Area', 'PixelList');
        % This is necessary if object is split by by cutting out local bbox
        if length(obj1Prop) > 1
            comWrtEc = sum(bsxfun(@times, bsxfun(@minus, cat(1,obj1Prop(:).Centroid), par.kdbBorder), par.settings.scale) .^ 2,2);
            [~, comIdx] = min(comWrtEc);
            obj1Prop = obj1Prop(comIdx);
        end
        obj2Prop = regionprops(obj2, 'Centroid', 'Area', 'PixelList');
        % This as well
        if length(obj2Prop) > 1
            comWrtEc = sum(bsxfun(@times, bsxfun(@minus, cat(1,obj2Prop(:).Centroid), par.kdbBorder), par.settings.scale) .^ 2,2);
            [~, comIdx] = min(comWrtEc);
            obj2Prop = obj2Prop(comIdx);
        end
        % Calculate mission direction in nm
        direction = (obj2Prop.Centroid([2 1 3]) - obj1Prop.Centroid([2 1 3])) .* par.settings.scale;
        direction = direction ./ norm(direction);
        % First: Find a orthogonal vector & 2 more vec for ONB (uniformly sampled from vectors in orthogonal plane)
        [randomOrthVec, orth1, orth2] = findRandomOrthogonal(direction);
        % Calculate extend of start object wrt to error center in nm 
        pixelPosWrtErrorCenter = bsxfun(@times, bsxfun(@minus, obj1Prop.PixelList(:,[2 1 3]), par.kdbBorder), par.settings.scale);
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
            v(j).missions(i).isControl = logical(par.controlFlag);
            v(j).missions(i).isTutorial = logical(par.tutorialFlag);
            % Start and end segment used for rendering (section-local)
            v(j).missions(i).start.id = edges(i,1);
            v(j).missions(i).end.id = edges(i,2);
            % error center in voxel (global)
            v(j).missions(i).errorCenter = errorCenter + par.bboxBig(:,1)';
            % Calculate CoM of objects in voxel
            v(j).missions(i).start.CoM = obj1Prop.Centroid([2 1 3]) + v(j).missions(i).errorCenter - par.kdbBorder;
            v(j).missions(i).end.CoM = obj2Prop.Centroid([2 1 3]) + v(j).missions(i).errorCenter - par.kdbBorder;
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
    save(par.counter, 'counter');

    % Write mission json (named according to input file)
    missionStruct.sectionId = 'section1';
    missionStruct.viewports = v;
    savejson('', missionStruct, 'FileName', [rootFolder 'missions.json'], 'ParseLogical', 1);

    % save everything to synced folder for problem inspection if flag is set
    if par.saveForProblemInspector
        raw = loadRawData(par.raw.root, par.raw.prefix, par.bboxBig, 0); 
        save(['/zdata/manuel/sync/problemInspector/' 'kDbResults' strrep(par.start, '/', '-') '.mat'], 'par', 'raw', 'seg', 'mito', 'v', '-v7.3');
    end

end

function [randomVec, orth1, orth2] = findRandomOrthogonal(v)
    % returns random orthogonal vector and and the 2nd and 3rd vector for ONB 
    v = v ./ norm(v);
    if all(v == [0 0 1])
        orth1 = [1 0 -1].*v([3 2 1]);  
    else
        orth1 = [1 -1 0].*v([2 1 3]);
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

